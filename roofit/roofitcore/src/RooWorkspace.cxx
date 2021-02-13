/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooWorkspace.cxx
\class RooWorkspace
\ingroup Roofitcore

The RooWorkspace is a persistable container for RooFit projects. A workspace
can contain and own variables, p.d.f.s, functions and datasets. All objects
that live in the workspace are owned by the workspace. The `import()` method
enforces consistency of objects upon insertion into the workspace (e.g. no
duplicate object with the same name are allowed) and makes sure all objects
in the workspace are connected to each other. Easy accessor methods like
`pdf()`, `var()` and `data()` allow to refer to the contents of the workspace by
object name. The entire RooWorkspace can be saved into a ROOT TFile and organises
the consistent streaming of its contents without duplication.
If a RooWorkspace contains custom classes, i.e. classes not in the
ROOT distribution, portability of workspaces can be enhanced by
storing the source code of those classes in the workspace as well.
This process is also organized by the workspace through the
`importClassCode()` method.

### Seemingly random crashes when reading large workspaces
When reading or loading workspaces with deeply nested PDFs, one can encounter
ouf-of-memory errors if the stack size is too small. This manifests in crashes
at seemingly random locations, or in the process silently ending.
Unfortunately, ROOT neither recover from this situation, nor warn or give useful
instructions. When suspecting to have run out of stack memory, check
```
ulimit -s
```
and try reading again.
**/

#include "RooWorkspace.h"
#include "RooHelpers.h"
#include "RooFit.h"
#include "RooWorkspaceHandle.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooAbsData.h"
#include "RooCmdConfig.h"
#include "RooMsgService.h"
#include "RooConstVar.h"
#include "RooResolutionModel.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "TBuffer.h"
#include "TInterpreter.h"
#include "TClassTable.h"
#include "TBaseClass.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "RooFactoryWSTool.h"
#include "RooAbsStudy.h"
#include "RooTObjWrap.h"
#include "RooAbsOptTestStatistic.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TClass.h"
#include "strlcpy.h"

#include <map>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

ClassImp(RooWorkspace);

////////////////////////////////////////////////////////////////////////////////

ClassImp(RooWorkspace::CodeRepo);

////////////////////////////////////////////////////////////////////////////////

ClassImp(RooWorkspace::WSDir);

list<string> RooWorkspace::_classDeclDirList ;
list<string> RooWorkspace::_classImplDirList ;
string RooWorkspace::_classFileExportDir = ".wscode.%s.%s" ;
Bool_t RooWorkspace::_autoClass = kFALSE ;


////////////////////////////////////////////////////////////////////////////////
/// Add `dir` to search path for class declaration (header) files. This is needed
/// to find class headers custom classes are imported into the workspace.
void RooWorkspace::addClassDeclImportDir(const char* dir)
{
  _classDeclDirList.push_back(dir) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Add `dir` to search path for class implementation (.cxx) files. This is needed
/// to find class headers custom classes are imported into the workspace.
void RooWorkspace::addClassImplImportDir(const char* dir)
{
  _classImplDirList.push_back(dir) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Specify the name of the directory in which embedded source
/// code is unpacked and compiled. The specified string may contain
/// one '%s' token which will be substituted by the workspace name

void RooWorkspace::setClassFileExportDir(const char* dir)
{
  if (dir) {
    _classFileExportDir = dir ;
  } else {
    _classFileExportDir = ".wscode.%s.%s" ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// If flag is true, source code of classes not the the ROOT distribution
/// is automatically imported if on object of such a class is imported
/// in the workspace

void RooWorkspace::autoImportClassCode(Bool_t flag)
{
  _autoClass = flag ;
}



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooWorkspace::RooWorkspace() : _classes(this), _dir(nullptr), _factory(nullptr), _doExport(kFALSE), _openTrans(kFALSE)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Construct empty workspace with given name and title

RooWorkspace::RooWorkspace(const char* name, const char* title) :
  TNamed(name,title?title:name), _classes(this), _dir(nullptr), _factory(nullptr), _doExport(kFALSE), _openTrans(kFALSE)
{
}


RooWorkspace::RooWorkspace(const char* name, Bool_t doCINTExport)  :
  TNamed(name,name), _classes(this), _dir(nullptr), _factory(nullptr), _doExport(kFALSE), _openTrans(kFALSE)
{
  // Construct empty workspace with given name and option to export reference to all workspace contents to a CINT namespace with the same name
  if (doCINTExport) {
    exportToCint(name) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Workspace copy constructor

RooWorkspace::RooWorkspace(const RooWorkspace& other) :
  TNamed(other), _uuid(other._uuid), _classes(other._classes,this), _dir(nullptr), _factory(nullptr), _doExport(kFALSE), _openTrans(kFALSE)
{
  // Copy owned nodes
  other._allOwnedNodes.snapshot(_allOwnedNodes,kTRUE) ;

  // Copy datasets
  TIterator* iter = other._dataList.MakeIterator() ;
  TObject* data2 ;
  while((data2=iter->Next())) {
    _dataList.Add(data2->Clone()) ;
  }
  delete iter ;

  // Copy snapshots
  TIterator* iter2 = other._snapshots.MakeIterator() ;
  RooArgSet* snap ;
  while((snap=(RooArgSet*)iter2->Next())) {
    RooArgSet* snapClone = (RooArgSet*) snap->snapshot() ;
    snapClone->setName(snap->GetName()) ;
    _snapshots.Add(snapClone) ;
  }
  delete iter2 ;

  // Copy named sets
  for (map<string,RooArgSet>::const_iterator iter3 = other._namedSets.begin() ; iter3 != other._namedSets.end() ; ++iter3) {
    // Make RooArgSet with equivalent content of this workspace
    RooArgSet* tmp = (RooArgSet*) _allOwnedNodes.selectCommon(iter3->second) ;
    _namedSets[iter3->first].add(*tmp) ;
    delete tmp ;
  }

  // Copy generic objects
  TIterator* iter4 = other._genObjects.MakeIterator() ;
  TObject* gobj ;
  while((gobj=iter4->Next())) {
    TObject *theClone = gobj->Clone();

    auto handle = dynamic_cast<RooWorkspaceHandle*>(theClone);
    if (handle) {
      handle->ReplaceWS(this);
    }

    _genObjects.Add(theClone);
  }
  delete iter4 ;

}



////////////////////////////////////////////////////////////////////////////////
/// Workspace destructor

RooWorkspace::~RooWorkspace()
{
  // Delete references to variables that were declared in CINT
  if (_doExport) {
    unExport() ;
  }

  // Delete contents
  _dataList.Delete() ;
  if (_dir) {
    delete _dir ;
  }
  _snapshots.Delete() ;

  // WVE named sets too?

  _genObjects.Delete() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Import a RooAbsArg or RooAbsData set from a workspace in a file. Filespec should be constructed as "filename:wspacename:objectname"
/// The arguments will be passed to the relevant import() or import(RooAbsData&, ...) import calls
/// \note From python, use `Import()`, since `import` is a reserved keyword.
Bool_t RooWorkspace::import(const char* fileSpec,
			    const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3,
			    const RooCmdArg& arg4, const RooCmdArg& arg5, const RooCmdArg& arg6,
			    const RooCmdArg& arg7, const RooCmdArg& arg8, const RooCmdArg& arg9)
{
  // Parse file/workspace/objectname specification
  std::vector<std::string> tokens = RooHelpers::tokenise(fileSpec, ":");

  // Check that parsing was successful
  if (tokens.size() != 3) {
    std::ostringstream stream;
    for (const auto& token : tokens) {
      stream << "\n\t" << token;
    }
    coutE(InputArguments) << "RooWorkspace(" << GetName() << ") ERROR in file specification, expecting 'filename:wsname:objname', but '" << fileSpec << "' given."
        << "\nTokens read are:" << stream.str() << endl;
    return kTRUE ;
  }

  const std::string& filename = tokens[0];
  const std::string& wsname = tokens[1];
  const std::string& objname = tokens[2];

  // Check that file can be opened
  TFile* f = TFile::Open(filename.c_str()) ;
  if (f==0) {
    coutE(InputArguments) << "RooWorkspace(" << GetName() << ") ERROR opening file " << filename << endl ;
    return 0 ;
  }

  // That that file contains workspace
  RooWorkspace* w = dynamic_cast<RooWorkspace*>(f->Get(wsname.c_str())) ;
  if (w==0) {
    coutE(InputArguments) << "RooWorkspace(" << GetName() << ") ERROR: No object named " << wsname << " in file " << filename
			  << " or object is not a RooWorkspace" << endl ;
    return 0 ;
  }

  // Check that workspace contains object and forward to appropriate import method
  RooAbsArg* warg = w->arg(objname.c_str()) ;
  if (warg) {
    Bool_t ret = import(*warg,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9) ;
    delete f ;
    return ret ;
  }
  RooAbsData* wdata = w->data(objname.c_str()) ;
  if (wdata) {
    Bool_t ret = import(*wdata,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9) ;
    delete f ;
    return ret ;
  }

  coutE(InputArguments) << "RooWorkspace(" << GetName() << ") ERROR: No RooAbsArg or RooAbsData object named " << objname
			<< " in workspace " << wsname << " in file " << filename << endl ;
  return kTRUE ;
}


////////////////////////////////////////////////////////////////////////////////
/// Import multiple RooAbsArg objects into workspace. For details on arguments see documentation
/// of import() method for single RooAbsArg
/// \note From python, use `Import()`, since `import` is a reserved keyword.
Bool_t RooWorkspace::import(const RooArgSet& args,
			    const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3,
			    const RooCmdArg& arg4, const RooCmdArg& arg5, const RooCmdArg& arg6,
			    const RooCmdArg& arg7, const RooCmdArg& arg8, const RooCmdArg& arg9)
{
  unique_ptr<TIterator> iter(args.createIterator()) ;
  RooAbsArg* oneArg ;
  Bool_t ret(kFALSE) ;
  while((oneArg=(RooAbsArg*)iter->Next())) {
    ret |= import(*oneArg,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9) ;
  }
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
///  Import a RooAbsArg object, e.g. function, p.d.f or variable into the workspace. This import function clones the input argument and will
///  own the clone. If a composite object is offered for import, e.g. a p.d.f with parameters and observables, the
///  complete tree of objects is imported. If any of the _variables_ of a composite object (parameters/observables) are already
///  in the workspace the imported p.d.f. is connected to the already existing variables. If any of the _function_ objects (p.d.f, formulas)
///  to be imported already exists in the workspace an error message is printed and the import of the entire tree of objects is cancelled.
///  Several optional arguments can be provided to modify the import procedure.
///
///  <table>
///  <tr><th> Accepted arguments
///  <tr><td> `RenameConflictNodes(const char* suffix)`   <td>  Add suffix to branch node name if name conflicts with existing node in workspace
///  <tr><td> `RenameAllNodes(const char* suffix)`    <td>  Add suffix to all branch node names including top level node.
///  <tr><td> `RenameAllVariables(const char* suffix)`    <td>  Add suffix to all variables of objects being imported.
///  <tr><td> `RenameAllVariablesExcept(const char* suffix, const char* exceptionList)`   <td>  Add suffix to all variables names, except ones listed
///  <tr><td> `RenameVariable(const char* inputName, const char* outputName)` <td>  Rename a single variable as specified upon import.
///  <tr><td> `RecycleConflictNodes()`    <td>  If any of the function objects to be imported already exist in the name space, connect the
///                            imported expression to the already existing nodes.
///                            \attention Use with care! If function definitions do not match, this alters the definition of your function upon import
///
///  <tr><td> `Silence()` <td>  Do not issue any info message
///  </table>
///
///  The RenameConflictNodes, RenameNodes and RecycleConflictNodes arguments are mutually exclusive. The RenameVariable argument can be repeated
///  as often as necessary to rename multiple variables. Alternatively, a single RenameVariable argument can be given with
///  two comma separated lists.
/// \note From python, use `Import()`, since `import` is a reserved keyword.
Bool_t RooWorkspace::import(const RooAbsArg& inArg,
    const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3,
    const RooCmdArg& arg4, const RooCmdArg& arg5, const RooCmdArg& arg6,
    const RooCmdArg& arg7, const RooCmdArg& arg8, const RooCmdArg& arg9)
{
  RooLinkedList args ;
  args.Add((TObject*)&arg1) ;
  args.Add((TObject*)&arg2) ;
  args.Add((TObject*)&arg3) ;
  args.Add((TObject*)&arg4) ;
  args.Add((TObject*)&arg5) ;
  args.Add((TObject*)&arg6) ;
  args.Add((TObject*)&arg7) ;
  args.Add((TObject*)&arg8) ;
  args.Add((TObject*)&arg9) ;

  // Select the pdf-specific commands
  RooCmdConfig pc(Form("RooWorkspace::import(%s)",GetName())) ;

  pc.defineString("conflictSuffix","RenameConflictNodes",0) ;
  pc.defineInt("renameConflictOrig","RenameConflictNodes",0,0) ;
  pc.defineString("allSuffix","RenameAllNodes",0) ;
  pc.defineString("allVarsSuffix","RenameAllVariables",0) ;
  pc.defineString("allVarsExcept","RenameAllVariables",1) ;
  pc.defineString("varChangeIn","RenameVar",0,"",kTRUE) ;
  pc.defineString("varChangeOut","RenameVar",1,"",kTRUE) ;
  pc.defineString("factoryTag","FactoryTag",0) ;
  pc.defineInt("useExistingNodes","RecycleConflictNodes",0,0) ;
  pc.defineInt("silence","Silence",0,0) ;
  pc.defineInt("noRecursion","NoRecursion",0,0) ;
  pc.defineMutex("RenameConflictNodes","RenameAllNodes") ;
  pc.defineMutex("RenameConflictNodes","RecycleConflictNodes") ;
  pc.defineMutex("RenameAllNodes","RecycleConflictNodes") ;
  pc.defineMutex("RenameVariable","RenameAllVariables") ;

  // Process and check varargs
  pc.process(args) ;
  if (!pc.ok(kTRUE)) {
    return kTRUE ;
  }

  // Decode renaming logic into suffix string and boolean for conflictOnly mode
  const char* suffixC = pc.getString("conflictSuffix") ;
  const char* suffixA = pc.getString("allSuffix") ;
  const char* suffixV = pc.getString("allVarsSuffix") ;
  const char* exceptVars = pc.getString("allVarsExcept") ;
  const char* varChangeIn = pc.getString("varChangeIn") ;
  const char* varChangeOut = pc.getString("varChangeOut") ;
  Bool_t renameConflictOrig = pc.getInt("renameConflictOrig") ;
  Int_t useExistingNodes = pc.getInt("useExistingNodes") ;
  Int_t silence = pc.getInt("silence") ;
  Int_t noRecursion = pc.getInt("noRecursion") ;


  // Turn zero length strings into null pointers
  if (suffixC && strlen(suffixC)==0) suffixC = 0 ;
  if (suffixA && strlen(suffixA)==0) suffixA = 0 ;

  Bool_t conflictOnly = suffixA ? kFALSE : kTRUE ;
  const char* suffix = suffixA ? suffixA : suffixC ;

  // Process any change in variable names
  map<string,string> varMap ;
  if (strlen(varChangeIn)>0) {

    // Parse comma separated lists into map<string,string>
    const std::vector<std::string> tokIn = RooHelpers::tokenise(varChangeIn, ", ");
    const std::vector<std::string> tokOut = RooHelpers::tokenise(varChangeOut, ", ");
    for (unsigned int i=0; i < tokIn.size(); ++i) {
      varMap.insert(std::make_pair(tokIn[i], tokOut[i]));
    }

    assert(tokIn.size() == tokOut.size());
  }

  // Process RenameAllVariables argument if specified
  // First convert exception list if provided
  std::set<string> exceptVarNames ;
  if (exceptVars && strlen(exceptVars)) {
    const std::vector<std::string> toks = RooHelpers::tokenise(exceptVars, ", ");
    exceptVarNames.insert(toks.begin(), toks.end());
  }

  if (suffixV != 0 && strlen(suffixV)>0) {
    RooArgSet* vars = inArg.getVariables() ;
    for (const auto v : *vars) {
      if (exceptVarNames.find(v->GetName())==exceptVarNames.end()) {
        varMap[v->GetName()] = Form("%s_%s",v->GetName(),suffixV) ;
      }
    }
    delete vars ;
  }

  // Scan for overlaps with current contents
  RooAbsArg* wsarg = _allOwnedNodes.find(inArg.GetName()) ;

  // Check for factory specification match
  const char* tagIn = inArg.getStringAttribute("factory_tag") ;
  const char* tagWs = wsarg ? wsarg->getStringAttribute("factory_tag") : 0 ;
  Bool_t factoryMatch = (tagIn && tagWs && !strcmp(tagIn,tagWs)) ;
  if (factoryMatch) {
    ((RooAbsArg&)inArg).setAttribute("RooWorkspace::Recycle") ;
  }

  if (!suffix && wsarg && !useExistingNodes && !(inArg.isFundamental() && varMap[inArg.GetName()]!="")) {
    if (!factoryMatch) {
      if (wsarg!=&inArg) {
        coutE(ObjectHandling) << "RooWorkSpace::import(" << GetName() << ") ERROR importing object named " << inArg.GetName()
			          << ": another instance with same name already in the workspace and no conflict resolution protocol specified" << endl ;
        return kTRUE ;
      } else {
        if (!silence) {
          coutI(ObjectHandling) << "RooWorkSpace::import(" << GetName() << ") Object " << inArg.GetName() << " is already in workspace!" << endl ;
        }
        return kTRUE ;
      }
    } else {
      coutI(ObjectHandling) << "RooWorkSpace::import(" << GetName() << ") Recycling existing object " << inArg.GetName() << " created with identical factory specification" << endl ;
    }
  }

  // Make list of conflicting nodes
  RooArgSet conflictNodes ;
  RooArgSet branchSet ;
  if (noRecursion) {
    branchSet.add(inArg) ;
  } else {
    inArg.branchNodeServerList(&branchSet) ;
  }

  for (const auto branch : branchSet) {
    RooAbsArg* wsbranch = _allOwnedNodes.find(branch->GetName()) ;
    if (wsbranch && wsbranch!=branch && !branch->getAttribute("RooWorkspace::Recycle") && !useExistingNodes) {
      conflictNodes.add(*branch) ;
    }
  }

  // Terminate here if there are conflicts and no resolution protocol
  if (conflictNodes.getSize()>0 && !suffix && !useExistingNodes) {
    coutE(ObjectHandling) << "RooWorkSpace::import(" << GetName() << ") ERROR object named " << inArg.GetName() << ": component(s) "
        << conflictNodes << " already in the workspace and no conflict resolution protocol specified" << endl ;
    return kTRUE ;
  }

  // Now create a working copy of the incoming object tree
  RooArgSet* cloneSet = new RooArgSet();
  cloneSet->useHashMapForFind(true); // Accelerate finding
  RooArgSet(inArg).snapshot(*cloneSet, !noRecursion);
  RooAbsArg* cloneTop = cloneSet->find(inArg.GetName()) ;

  // Mark all nodes for renaming if we are not in conflictOnly mode
  if (!conflictOnly) {
    conflictNodes.removeAll() ;
    conflictNodes.add(branchSet) ;
  }

  // Mark nodes that are to be renamed with special attribute
  string topName2 = cloneTop->GetName() ;
  if (!renameConflictOrig) {
    // Mark all nodes to be imported for renaming following conflict resolution protocol
    for (const auto cnode : conflictNodes) {
      RooAbsArg* cnode2 = cloneSet->find(cnode->GetName()) ;
      string origName = cnode2->GetName() ;
      cnode2->SetName(Form("%s_%s",cnode2->GetName(),suffix)) ;
      cnode2->SetTitle(Form("%s (%s)",cnode2->GetTitle(),suffix)) ;
      string tag = Form("ORIGNAME:%s",origName.c_str()) ;
      cnode2->setAttribute(tag.c_str()) ;
      if (!cnode2->getStringAttribute("origName")) {
        string tag2 = Form("%s",origName.c_str()) ;
        cnode2->setStringAttribute("origName",tag2.c_str()) ;
      }

      // Save name of new top level node for later use
      if (cnode2==cloneTop) {
        topName2 = cnode2->GetName() ;
      }

      if (!silence) {
        coutI(ObjectHandling) << "RooWorkspace::import(" << GetName()
			          << ") Resolving name conflict in workspace by changing name of imported node  "
			          << origName << " to " << cnode2->GetName() << endl ;
      }
    }
  } else {

    // Rename all nodes already in the workspace to 'clear the way' for the imported nodes
    for (const auto cnode : conflictNodes) {

      string origName = cnode->GetName() ;
      RooAbsArg* wsnode = _allOwnedNodes.find(origName.c_str()) ;
      if (wsnode) {

        if (!wsnode->getStringAttribute("origName")) {
          wsnode->setStringAttribute("origName",wsnode->GetName()) ;
        }

        if (!_allOwnedNodes.find(Form("%s_%s",cnode->GetName(),suffix))) {
          wsnode->SetName(Form("%s_%s",cnode->GetName(),suffix)) ;
          wsnode->SetTitle(Form("%s (%s)",cnode->GetTitle(),suffix)) ;
        } else {
          // Name with suffix already taken, add additional suffix
          for (unsigned int n=1; true; ++n) {
            string newname = Form("%s_%s_%d",cnode->GetName(),suffix,n) ;
            if (!_allOwnedNodes.find(newname.c_str())) {
              wsnode->SetName(newname.c_str()) ;
              wsnode->SetTitle(Form("%s (%s %d)",cnode->GetTitle(),suffix,n)) ;
              break ;
            }
          }
        }
        if (!silence) {
          coutI(ObjectHandling) << "RooWorkspace::import(" << GetName()
				    << ") Resolving name conflict in workspace by changing name of original node "
				    << origName << " to " << wsnode->GetName() << endl ;
        }
      } else {
        coutW(ObjectHandling) << "RooWorkspce::import(" << GetName() << ") Internal error: expected to find existing node "
            << origName << " to be renamed, but didn't find it..." << endl ;
      }

    }
  }

  // Process any change in variable names
  if (strlen(varChangeIn)>0 || (suffixV && strlen(suffixV)>0)) {

    // Process all changes in variable names
    for (const auto cnode : *cloneSet) {

      if (varMap.find(cnode->GetName())!=varMap.end()) {
        string origName = cnode->GetName() ;
        cnode->SetName(varMap[cnode->GetName()].c_str()) ;
        string tag = Form("ORIGNAME:%s",origName.c_str()) ;
        cnode->setAttribute(tag.c_str()) ;
        if (!cnode->getStringAttribute("origName")) {
          string tag2 = Form("%s",origName.c_str()) ;
          cnode->setStringAttribute("origName",tag2.c_str()) ;
        }

        if (!silence) {
          coutI(ObjectHandling) << "RooWorkspace::import(" << GetName() << ") Changing name of variable "
              << origName << " to " << cnode->GetName() << " on request" << endl ;
        }

        if (cnode==cloneTop) {
          topName2 = cnode->GetName() ;
        }

      }
    }
  }

  // Now clone again with renaming effective
  RooArgSet* cloneSet2 = new RooArgSet();
  cloneSet2->useHashMapForFind(true); // Faster finding
  RooArgSet(*cloneTop).snapshot(*cloneSet2, !noRecursion);
  RooAbsArg* cloneTop2 = cloneSet2->find(topName2.c_str()) ;

  // Make final check list of conflicting nodes
  RooArgSet conflictNodes2 ;
  RooArgSet branchSet2 ;
  for (const auto branch2 : branchSet2) {
    if (_allOwnedNodes.find(branch2->GetName())) {
      conflictNodes2.add(*branch2) ;
    }
  }

  // Terminate here if there are conflicts and no resolution protocol
  if (conflictNodes2.getSize()) {
    coutE(ObjectHandling) << "RooWorkSpace::import(" << GetName() << ") ERROR object named " << inArg.GetName() << ": component(s) "
        << conflictNodes2 << " cause naming conflict after conflict resolution protocol was executed" << endl ;
    return kTRUE ;
  }

  // Perform any auxiliary imports at this point
  for (const auto node : *cloneSet2) {
    if (node->importWorkspaceHook(*this)) {
      coutE(ObjectHandling) << "RooWorkSpace::import(" << GetName() << ") ERROR object named " << node->GetName()
			        << " has an error in importing in one or more of its auxiliary objects, aborting" << endl ;
      return kTRUE ;
    }
  }

  RooArgSet recycledNodes ;
  RooArgSet nodesToBeDeleted ;
  for (const auto node : *cloneSet2) {
    if (_autoClass) {
      if (!_classes.autoImportClass(node->IsA())) {
        coutW(ObjectHandling) << "RooWorkspace::import(" << GetName() << ") WARNING: problems import class code of object "
            << node->IsA()->GetName() << "::" << node->GetName() << ", reading of workspace will require external definition of class" << endl ;
      }
    }

    // Point expensiveObjectCache to copy in this workspace
    RooExpensiveObjectCache& oldCache = node->expensiveObjectCache() ;
    node->setExpensiveObjectCache(_eocache) ;
    _eocache.importCacheObjects(oldCache,node->GetName(),kTRUE) ;

    // Check if node is already in workspace (can only happen for variables or identical instances, unless RecycleConflictNodes is specified)
    RooAbsArg* wsnode = _allOwnedNodes.find(node->GetName()) ;

    if (wsnode) {
      // Do not import node, add not to list of nodes that require reconnection
      if (!silence && useExistingNodes) {
        coutI(ObjectHandling) << "RooWorkspace::import(" << GetName() << ") using existing copy of " << node->IsA()->GetName()
			          << "::" << node->GetName() << " for import of " << cloneTop2->IsA()->GetName() << "::"
			          << cloneTop2->GetName() << endl ;
      }
      recycledNodes.add(*_allOwnedNodes.find(node->GetName())) ;

      // Delete clone of incoming node
      nodesToBeDeleted.addOwned(*node) ;

      //cout << "WV: recycling existing node " << existingNode << " = " << existingNode->GetName() << " for imported node " << node << endl ;

    } else {
      // Import node
      if (!silence) {
        coutI(ObjectHandling) << "RooWorkspace::import(" << GetName() << ") importing " << node->IsA()->GetName() << "::"
            << node->GetName() << endl ;
      }
      _allOwnedNodes.addOwned(*node) ;
      if (_openTrans) {
        _sandboxNodes.add(*node) ;
      } else {
        if (_dir && node->IsA() != RooConstVar::Class()) {
          _dir->InternalAppend(node) ;
        }
        if (_doExport && node->IsA() != RooConstVar::Class()) {
          exportObj(node) ;
        }
      }
    }
  }

  // Release working copy
  // no need to do a safe list erase since it was generated from a snapshot
  // just take ownership and delte elements by hand
  cloneSet->releaseOwnership();
  for (auto cloneNode : *cloneSet){
    delete cloneNode;
  }
  delete cloneSet ;

  // Reconnect any nodes that need to be
  if (recycledNodes.getSize()>0) {
    for (const auto node : *cloneSet2) {
      node->redirectServers(recycledNodes) ;
    }
  }

  cloneSet2->releaseOwnership() ;
  delete cloneSet2 ;

  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
///  Import a dataset (RooDataSet or RooDataHist) into the work space. The workspace will contain a copy of the data.
///  The dataset and its variables can be renamed upon insertion with the options below
///
///  <table>
/// <tr><th> Accepted arguments
/// <tr><td> `Rename(const char* suffix)` <td> Rename dataset upon insertion
/// <tr><td> `RenameVariable(const char* inputName, const char* outputName)` <td> Change names of observables in dataset upon insertion
/// <tr><td> `Silence` <td> Be quiet, except in case of errors
/// \note From python, use `Import()`, since `import` is a reserved keyword.
Bool_t RooWorkspace::import(RooAbsData& inData,
			    const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3,
			    const RooCmdArg& arg4, const RooCmdArg& arg5, const RooCmdArg& arg6,
			    const RooCmdArg& arg7, const RooCmdArg& arg8, const RooCmdArg& arg9)

{

  RooLinkedList args ;
  args.Add((TObject*)&arg1) ;
  args.Add((TObject*)&arg2) ;
  args.Add((TObject*)&arg3) ;
  args.Add((TObject*)&arg4) ;
  args.Add((TObject*)&arg5) ;
  args.Add((TObject*)&arg6) ;
  args.Add((TObject*)&arg7) ;
  args.Add((TObject*)&arg8) ;
  args.Add((TObject*)&arg9) ;

  // Select the pdf-specific commands
  RooCmdConfig pc(Form("RooWorkspace::import(%s)",GetName())) ;

  pc.defineString("dsetName","Rename",0,"") ;
  pc.defineString("varChangeIn","RenameVar",0,"",kTRUE) ;
  pc.defineString("varChangeOut","RenameVar",1,"",kTRUE) ;
  pc.defineInt("embedded","Embedded",0,0) ;
  pc.defineInt("silence","Silence",0,0) ;

  // Process and check varargs
  pc.process(args) ;
  if (!pc.ok(kTRUE)) {
    return kTRUE ;
  }

  // Decode renaming logic into suffix string and boolean for conflictOnly mode
  const char* dsetName = pc.getString("dsetName") ;
  const char* varChangeIn = pc.getString("varChangeIn") ;
  const char* varChangeOut = pc.getString("varChangeOut") ;
  Bool_t embedded = pc.getInt("embedded") ;
  Int_t silence = pc.getInt("silence") ;

  if (!silence)
    coutI(ObjectHandling) << "RooWorkspace::import(" << GetName() << ") importing dataset " << inData.GetName() << endl ;

  // Transform emtpy string into null pointer
  if (dsetName && strlen(dsetName)==0) {
    dsetName=0 ;
  }

  RooLinkedList& dataList = embedded ? _embeddedDataList : _dataList ;
  if (dataList.GetSize() > 50 && dataList.getHashTableSize() == 0) {
    // When the workspaces get larger, traversing the linked list becomes a bottleneck:
    dataList.setHashTableSize(200);
  }

  // Check that no dataset with target name already exists
  if (dsetName && dataList.FindObject(dsetName)) {
    coutE(ObjectHandling) << "RooWorkspace::import(" << GetName() << ") ERROR dataset with name " << dsetName << " already exists in workspace, import aborted" << endl ;
    return kTRUE ;
  }
  if (!dsetName && dataList.FindObject(inData.GetName())) {
    coutE(ObjectHandling) << "RooWorkspace::import(" << GetName() << ") ERROR dataset with name " << inData.GetName() << " already exists in workspace, import aborted" << endl ;
    return kTRUE ;
  }

  // Rename dataset if required
  RooAbsData* clone ;
  if (dsetName) {
    if (!silence)
      coutI(ObjectHandling) << "RooWorkSpace::import(" << GetName() << ") changing name of dataset from  " << inData.GetName() << " to " << dsetName << endl ;
    clone = (RooAbsData*) inData.Clone(dsetName) ;
  } else {
    clone = (RooAbsData*) inData.Clone(inData.GetName()) ;
  }


  // Process any change in variable names
  if (strlen(varChangeIn)>0) {
    // Parse comma separated lists of variable name changes
    const std::vector<std::string> tokIn  = RooHelpers::tokenise(varChangeIn, ",");
    const std::vector<std::string> tokOut = RooHelpers::tokenise(varChangeOut, ",");
    for (unsigned int i=0; i < tokIn.size(); ++i) {
      if (!silence)
        coutI(ObjectHandling) << "RooWorkSpace::import(" << GetName() << ") changing name of dataset observable " << tokIn[i] << " to " << tokOut[i] << endl ;
      clone->changeObservableName(tokIn[i].c_str(), tokOut[i].c_str());
    }
  }

  // Now import the dataset observables, unless dataset is embedded
  RooAbsArg* carg ;
  if (!embedded) {
    TIterator* iter = clone->get()->createIterator() ;
    while((carg=(RooAbsArg*)iter->Next())) {
      if (!arg(carg->GetName())) {
	import(*carg) ;
      }
    }
    delete iter ;
  }

  dataList.Add(clone) ;
  if (_dir) {
    _dir->InternalAppend(clone) ;
  }
  if (_doExport) {
    exportObj(clone) ;
  }

  // Set expensive object cache of dataset internal buffers to that of workspace
  RooFIter iter2 = clone->get()->fwdIterator() ;
  while ((carg=iter2.next())) {
    carg->setExpensiveObjectCache(expensiveObjectCache()) ;
  }


  return kFALSE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Define a named RooArgSet with given constituents. If importMissing is true, any constituents
/// of aset that are not in the workspace will be imported, otherwise an error is returned
/// for missing components

Bool_t RooWorkspace::defineSet(const char* name, const RooArgSet& aset, Bool_t importMissing)
{
  // Check if set was previously defined, if so print warning
  map<string,RooArgSet>::iterator i = _namedSets.find(name) ;
  if (i!=_namedSets.end()) {
    coutW(InputArguments) << "RooWorkspace::defineSet(" << GetName() << ") WARNING redefining previously defined named set " << name << endl ;
  }

  RooArgSet wsargs ;

  // Check all constituents of provided set
  TIterator* iter = aset.createIterator() ;
  RooAbsArg* sarg ;
  while((sarg=(RooAbsArg*)iter->Next())) {
    // If missing, either import or report error
    if (!arg(sarg->GetName())) {
      if (importMissing) {
	import(*sarg) ;
      } else {
	coutE(InputArguments) << "RooWorkspace::defineSet(" << GetName() << ") ERROR set constituent \"" << sarg->GetName()
			      << "\" is not in workspace and importMissing option is disabled" << endl ;
	return kTRUE ;
      }
    }
    wsargs.add(*arg(sarg->GetName())) ;
  }
  delete iter ;

  // Install named set
  _namedSets[name].removeAll() ;
  _namedSets[name].add(wsargs) ;

  return kFALSE ;
}

//_____________________________________________________________________________
Bool_t RooWorkspace::defineSetInternal(const char *name, const RooArgSet &aset)
{
   // Define a named RooArgSet with given constituents. If importMissing is true, any constituents
   // of aset that are not in the workspace will be imported, otherwise an error is returned
   // for missing components

   // Check if set was previously defined, if so print warning
   map<string, RooArgSet>::iterator i = _namedSets.find(name);
   if (i != _namedSets.end()) {
      coutW(InputArguments) << "RooWorkspace::defineSet(" << GetName()
                            << ") WARNING redefining previously defined named set " << name << endl;
   }

   // Install named set
   _namedSets[name].removeAll();
   _namedSets[name].add(aset);

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Define a named set in the work space through a comma separated list of
/// names of objects already in the workspace

Bool_t RooWorkspace::defineSet(const char* name, const char* contentList)
{
  // Check if set was previously defined, if so print warning
  map<string,RooArgSet>::iterator i = _namedSets.find(name) ;
  if (i!=_namedSets.end()) {
    coutW(InputArguments) << "RooWorkspace::defineSet(" << GetName() << ") WARNING redefining previously defined named set " << name << endl ;
  }

  RooArgSet wsargs ;

  // Check all constituents of provided set
  for (const std::string& token : RooHelpers::tokenise(contentList, ",")) {
    // If missing, either import or report error
    if (!arg(token.c_str())) {
      coutE(InputArguments) << "RooWorkspace::defineSet(" << GetName() << ") ERROR proposed set constituent \"" << token
			    << "\" is not in workspace" << endl ;
      return kTRUE ;
    }
    wsargs.add(*arg(token.c_str())) ;
  }

  // Install named set
  _namedSets[name].removeAll() ;
  _namedSets[name].add(wsargs) ;

  return kFALSE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Define a named set in the work space through a comma separated list of
/// names of objects already in the workspace

Bool_t RooWorkspace::extendSet(const char* name, const char* newContents)
{
  RooArgSet wsargs ;

  // Check all constituents of provided set
  for (const std::string& token : RooHelpers::tokenise(newContents, ",")) {
    // If missing, either import or report error
    if (!arg(token.c_str())) {
      coutE(InputArguments) << "RooWorkspace::defineSet(" << GetName() << ") ERROR proposed set constituent \"" << token
			    << "\" is not in workspace" << endl ;
      return kTRUE ;
    }
    wsargs.add(*arg(token.c_str())) ;
  }

  // Extend named set
  _namedSets[name].add(wsargs,kTRUE) ;

  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return pointer to previously defined named set with given nmame
/// If no such set is found a null pointer is returned

const RooArgSet* RooWorkspace::set(const char* name)
{
  map<string,RooArgSet>::iterator i = _namedSets.find(name) ;
  return (i!=_namedSets.end()) ? &(i->second) : 0 ;
}




////////////////////////////////////////////////////////////////////////////////
/// Rename set to a new name

Bool_t RooWorkspace::renameSet(const char* name, const char* newName)
{
  // First check if set exists
  if (!set(name)) {
    coutE(InputArguments) << "RooWorkspace::renameSet(" << GetName() << ") ERROR a set with name " << name
			  << " does not exist" << endl ;
    return kTRUE ;
  }

  // Check if no set exists with new name
  if (set(newName)) {
    coutE(InputArguments) << "RooWorkspace::renameSet(" << GetName() << ") ERROR a set with name " << newName
			  << " already exists" << endl ;
    return kTRUE ;
  }

  // Copy entry under 'name' to 'newName'
  _namedSets[newName].add(_namedSets[name]) ;

  // Remove entry under old name
  _namedSets.erase(name) ;

  return kFALSE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Remove a named set from the workspace

Bool_t RooWorkspace::removeSet(const char* name)
{
  // First check if set exists
  if (!set(name)) {
    coutE(InputArguments) << "RooWorkspace::removeSet(" << GetName() << ") ERROR a set with name " << name
			  << " does not exist" << endl ;
    return kTRUE ;
  }

  // Remove set with given name
  _namedSets.erase(name) ;

  return kFALSE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Open an import transaction operations. Returns kTRUE if successful, kFALSE
/// if there is already an ongoing transaction

Bool_t RooWorkspace::startTransaction()
{
  // Check that there was no ongoing transaction
  if (_openTrans) {
    return kFALSE ;
  }

  // Open transaction
  _openTrans = kTRUE ;
  return kTRUE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Cancel an ongoing import transaction. All objects imported since startTransaction()
/// will be removed and the transaction will be terminated. Return kTRUE if cancel operation
/// succeeds, return kFALSE if there was no open transaction

Bool_t RooWorkspace::cancelTransaction()
{
  // Check that there is an ongoing transaction
  if (!_openTrans) {
    return kFALSE ;
  }

  // Delete all objects in the sandbox
  TIterator* iter = _sandboxNodes.createIterator() ;
  RooAbsArg* tmpArg ;
  while((tmpArg=(RooAbsArg*)iter->Next())) {
    _allOwnedNodes.remove(*tmpArg) ;
  }
  delete iter ;
  _sandboxNodes.removeAll() ;

  // Mark transaction as finished
  _openTrans = kFALSE ;

  return kTRUE ;
}

Bool_t RooWorkspace::commitTransaction()
{
  // Commit an ongoing import transaction. Returns kTRUE if commit succeeded,
  // return kFALSE if there was no ongoing transaction

  // Check that there is an ongoing transaction
  if (!_openTrans) {
    return kFALSE ;
  }

  // Publish sandbox nodes in directory and/or CINT if requested
  TIterator* iter = _sandboxNodes.createIterator() ;
  RooAbsArg* sarg ;
  while((sarg=(RooAbsArg*)iter->Next())) {
    if (_dir && sarg->IsA() != RooConstVar::Class()) {
      _dir->InternalAppend(sarg) ;
    }
    if (_doExport && sarg->IsA() != RooConstVar::Class()) {
      exportObj(sarg) ;
    }
  }
  delete iter ;

  // Remove all committed objects from the sandbox
  _sandboxNodes.removeAll() ;

  // Mark transaction as finished
  _openTrans = kFALSE ;

  return kTRUE ;
}




////////////////////////////////////////////////////////////////////////////////

Bool_t RooWorkspace::importClassCode(TClass* theClass, Bool_t doReplace)
{
  return _classes.autoImportClass(theClass,doReplace) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Inport code of all classes in the workspace that have a class name
/// that matches pattern 'pat' and which are not found to be part of
/// the standard ROOT distribution. If doReplace is true any existing
/// class code saved in the workspace is replaced

Bool_t RooWorkspace::importClassCode(const char* pat, Bool_t doReplace)
{
  Bool_t ret(kTRUE) ;

  TRegexp re(pat,kTRUE) ;
  TIterator* iter = componentIterator() ;
  RooAbsArg* carg ;
  while((carg=(RooAbsArg*)iter->Next())) {
    TString className = carg->IsA()->GetName() ;
    if (className.Index(re)>=0 && !_classes.autoImportClass(carg->IsA(),doReplace)) {
      coutW(ObjectHandling) << "RooWorkspace::import(" << GetName() << ") WARNING: problems import class code of object "
			    << carg->IsA()->GetName() << "::" << carg->GetName() << ", reading of workspace will require external definition of class" << endl ;
      ret = kFALSE ;
    }
  }
  delete iter ;

  return ret ;
}





////////////////////////////////////////////////////////////////////////////////
/// Save snapshot of values and attributes (including "Constant") of given parameters.
/// \param[in] name Name of the snapshot.
/// \param[in] paramNames Comma-separated list of parameter names to be snapshot.
Bool_t RooWorkspace::saveSnapshot(const char* name, const char* paramNames)
{
  return saveSnapshot(name,argSet(paramNames),kFALSE) ;
}





////////////////////////////////////////////////////////////////////////////////
/// Save snapshot of values and attributes (including "Constant") of parameters 'params'.
/// If importValues is FALSE, the present values from the object in the workspace are
/// saved. If importValues is TRUE, the values of the objects passed in the 'params'
/// argument are saved

Bool_t RooWorkspace::saveSnapshot(const char* name, const RooArgSet& params, Bool_t importValues)
{
  RooArgSet* actualParams = (RooArgSet*) _allOwnedNodes.selectCommon(params) ;
  RooArgSet* snapshot = (RooArgSet*) actualParams->snapshot() ;
  delete actualParams ;

  snapshot->setName(name) ;

  if (importValues) {
    *snapshot = params ;
  }

  RooArgSet* oldSnap = (RooArgSet*) _snapshots.FindObject(name) ;
  if (oldSnap) {
    coutI(ObjectHandling) << "RooWorkspace::saveSnaphot(" << GetName() << ") replacing previous snapshot with name " << name << endl ;
    _snapshots.Remove(oldSnap) ;
    delete oldSnap ;
  }

  _snapshots.Add(snapshot) ;

  return kTRUE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Load the values and attributes of the parameters in the snapshot saved with
/// the given name

Bool_t RooWorkspace::loadSnapshot(const char* name)
{
  RooArgSet* snap = (RooArgSet*) _snapshots.find(name) ;
  if (!snap) {
    coutE(ObjectHandling) << "RooWorkspace::loadSnapshot(" << GetName() << ") no snapshot with name " << name << " is available" << endl ;
    return kFALSE ;
  }

  RooArgSet* actualParams = (RooArgSet*) _allOwnedNodes.selectCommon(*snap) ;
  *actualParams = *snap ;
  delete actualParams ;

  return kTRUE ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the RooArgSet containing a snapshot of variables contained in the workspace
///
/// Note that the variables of the objects in the snapshots are **copies** of the
/// variables in the workspace. To load the values of a snapshot in the workspace
/// variables, use loadSnapshot() instead.

const RooArgSet* RooWorkspace::getSnapshot(const char* name) const
{
  RooArgSet* snap = (RooArgSet*) _snapshots.find(name) ;
  if (!snap) {
    coutE(ObjectHandling) << "RooWorkspace::loadSnapshot(" << GetName() << ") no snapshot with name " << name << " is available" << endl ;
    return 0 ;
  }

  return snap ;
}



// //_____________________________________________________________________________
// RooAbsPdf* RooWorkspace::joinPdf(const char* jointPdfName, const char* indexName, const char* inputMapping)
// {
//   // Join given list of p.d.f.s into a simultaneous p.d.f with given name. If the named index category
//   // does not exist, it is created.
//   //
//   //  Example : joinPdf("simPdf","expIndex","A=pdfA,B=pdfB") ;
//   //
//   //            will return a RooSimultaneous named 'simPdf' with index category 'expIndex' with
//   //            state names A and B. Pdf 'pdfA' will be associated with state A and pdf 'pdfB'
//   //            will be associated with state B
//   //
//   return 0 ;
// }

// //_____________________________________________________________________________
// RooAbsData* RooWorkspace::joinData(const char* jointDataName, const char* indexName, const char* inputMapping)
// {
//   // Join given list of dataset into a joint dataset for use with a simultaneous pdf
//   // (as e.g. created by joingPdf"
//   //
//   //  Example : joinData("simData","expIndex","A=dataA,B=dataB") ;
//   //
//   //            will return a RooDataSet named 'simData' that consist of the entries of both
//   //            dataA and dataB. An extra category column 'expIndex' is added that labels
//   //            each entry with state 'A' and 'B' to indicate the originating dataset
//   return 0 ;
// }


////////////////////////////////////////////////////////////////////////////////
/// Retrieve p.d.f (RooAbsPdf) with given name. A null pointer is returned if not found

RooAbsPdf* RooWorkspace::pdf(const char* name) const
{
  return dynamic_cast<RooAbsPdf*>(_allOwnedNodes.find(name)) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Retrieve function (RooAbsReal) with given name. Note that all RooAbsPdfs are also RooAbsReals. A null pointer is returned if not found.

RooAbsReal* RooWorkspace::function(const char* name) const
{
  return dynamic_cast<RooAbsReal*>(_allOwnedNodes.find(name)) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Retrieve real-valued variable (RooRealVar) with given name. A null pointer is returned if not found

RooRealVar* RooWorkspace::var(const char* name) const
{
  return dynamic_cast<RooRealVar*>(_allOwnedNodes.find(name)) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Retrieve discrete variable (RooCategory) with given name. A null pointer is returned if not found

RooCategory* RooWorkspace::cat(const char* name) const
{
  return dynamic_cast<RooCategory*>(_allOwnedNodes.find(name)) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Retrieve discrete function (RooAbsCategory) with given name. A null pointer is returned if not found

RooAbsCategory* RooWorkspace::catfunc(const char* name) const
{
  return dynamic_cast<RooAbsCategory*>(_allOwnedNodes.find(name)) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return RooAbsArg with given name. A null pointer is returned if none is found.

RooAbsArg* RooWorkspace::arg(const char* name) const
{
  return _allOwnedNodes.find(name) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return set of RooAbsArgs matching to given list of names

RooArgSet RooWorkspace::argSet(const char* nameList) const
{
  RooArgSet ret ;

  for (const std::string& token : RooHelpers::tokenise(nameList, ",")) {
    RooAbsArg* oneArg = arg(token.c_str()) ;
    if (oneArg) {
      ret.add(*oneArg) ;
    } else {
      coutE(InputArguments) << " RooWorkspace::argSet(" << GetName() << ") no RooAbsArg named \"" << token << "\" in workspace" << endl ;
    }
  }
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return fundamental (i.e. non-derived) RooAbsArg with given name. Fundamental types
/// are e.g. RooRealVar, RooCategory. A null pointer is returned if none is found.

RooAbsArg* RooWorkspace::fundArg(const char* name) const
{
  RooAbsArg* tmp = arg(name) ;
  if (!tmp) {
    return 0 ;
  }
  return tmp->isFundamental() ? tmp : 0 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Retrieve dataset (binned or unbinned) with given name. A null pointer is returned if not found

RooAbsData* RooWorkspace::data(const char* name) const
{
  return (RooAbsData*)_dataList.FindObject(name) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Retrieve dataset (binned or unbinned) with given name. A null pointer is returned if not found

RooAbsData* RooWorkspace::embeddedData(const char* name) const
{
  return (RooAbsData*)_embeddedDataList.FindObject(name) ;
}




////////////////////////////////////////////////////////////////////////////////
/// Return set with all variable objects

RooArgSet RooWorkspace::allVars() const
{
  RooArgSet ret ;

  // Split list of components in pdfs, functions and variables
  TIterator* iter = _allOwnedNodes.createIterator() ;
  RooAbsArg* parg ;
  while((parg=(RooAbsArg*)iter->Next())) {
    if (parg->IsA()->InheritsFrom(RooRealVar::Class())) {
      ret.add(*parg) ;
    }
  }
  delete iter ;

  return ret ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return set with all category objects

RooArgSet RooWorkspace::allCats() const
{
  RooArgSet ret ;

  // Split list of components in pdfs, functions and variables
  TIterator* iter = _allOwnedNodes.createIterator() ;
  RooAbsArg* parg ;
  while((parg=(RooAbsArg*)iter->Next())) {
    if (parg->IsA()->InheritsFrom(RooCategory::Class())) {
      ret.add(*parg) ;
    }
  }
  delete iter ;

  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return set with all function objects

RooArgSet RooWorkspace::allFunctions() const
{
  RooArgSet ret ;

  // Split list of components in pdfs, functions and variables
  TIterator* iter = _allOwnedNodes.createIterator() ;
  RooAbsArg* parg ;
  while((parg=(RooAbsArg*)iter->Next())) {
    if (parg->IsA()->InheritsFrom(RooAbsReal::Class()) &&
	!parg->IsA()->InheritsFrom(RooAbsPdf::Class()) &&
	!parg->IsA()->InheritsFrom(RooConstVar::Class()) &&
	!parg->IsA()->InheritsFrom(RooRealVar::Class())) {
      ret.add(*parg) ;
    }
  }

  return ret ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return set with all category function objects

RooArgSet RooWorkspace::allCatFunctions() const
{
  RooArgSet ret ;

  // Split list of components in pdfs, functions and variables
  TIterator* iter = _allOwnedNodes.createIterator() ;
  RooAbsArg* parg ;
  while((parg=(RooAbsArg*)iter->Next())) {
    if (parg->IsA()->InheritsFrom(RooAbsCategory::Class()) &&
	!parg->IsA()->InheritsFrom(RooCategory::Class())) {
      ret.add(*parg) ;
    }
  }
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return set with all resolution model objects

RooArgSet RooWorkspace::allResolutionModels() const
{
  RooArgSet ret ;

  // Split list of components in pdfs, functions and variables
  TIterator* iter = _allOwnedNodes.createIterator() ;
  RooAbsArg* parg ;
  while((parg=(RooAbsArg*)iter->Next())) {
    if (parg->IsA()->InheritsFrom(RooResolutionModel::Class())) {
      if (!((RooResolutionModel*)parg)->isConvolved()) {
	ret.add(*parg) ;
      }
    }
  }
  return ret ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return set with all probability density function objects

RooArgSet RooWorkspace::allPdfs() const
{
  RooArgSet ret ;

  // Split list of components in pdfs, functions and variables
  TIterator* iter = _allOwnedNodes.createIterator() ;
  RooAbsArg* parg ;
  while((parg=(RooAbsArg*)iter->Next())) {
    if (parg->IsA()->InheritsFrom(RooAbsPdf::Class()) &&
	!parg->IsA()->InheritsFrom(RooResolutionModel::Class())) {
      ret.add(*parg) ;
    }
  }
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return list of all dataset in the workspace

list<RooAbsData*> RooWorkspace::allData() const
{
  list<RooAbsData*> ret ;
  TIterator* iter = _dataList.MakeIterator() ;
  RooAbsData* dat ;
  while((dat=(RooAbsData*)iter->Next())) {
    ret.push_back(dat) ;
  }
  delete iter ;
  return ret ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return list of all dataset in the workspace

list<RooAbsData*> RooWorkspace::allEmbeddedData() const
{
  list<RooAbsData*> ret ;
  TIterator* iter = _embeddedDataList.MakeIterator() ;
  RooAbsData* dat ;
  while((dat=(RooAbsData*)iter->Next())) {
    ret.push_back(dat) ;
  }
  delete iter ;
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return list of all generic objects in the workspace

list<TObject*> RooWorkspace::allGenericObjects() const
{
  list<TObject*> ret ;
  TIterator* iter = _genObjects.MakeIterator() ;
  TObject* gobj ;
  while((gobj=(RooAbsData*)iter->Next())) {

    // If found object is wrapper, return payload
    if (gobj->IsA()==RooTObjWrap::Class()) {
      ret.push_back(((RooTObjWrap*)gobj)->obj()) ;
    } else {
      ret.push_back(gobj) ;
    }
  }
  delete iter ;
  return ret ;
}




////////////////////////////////////////////////////////////////////////////////
/// Import code of class 'tc' into the repository. If code is already in repository it is only imported
/// again if doReplace is false. The names and location of the source files is determined from the information
/// in TClass. If no location is found in the TClass information, the files are searched in the workspace
/// search path, defined by addClassDeclImportDir() and addClassImplImportDir() for declaration and implementation
/// files respectively. If files cannot be found, abort with error status, otherwise update the internal
/// class-to-file map and import the contents of the files, if they are not imported yet.

Bool_t RooWorkspace::CodeRepo::autoImportClass(TClass* tc, Bool_t doReplace)
{

  oocxcoutD(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo(" << _wspace->GetName() << ") request to import code of class " << tc->GetName() << endl ;

  // *** PHASE 1 *** Check if file needs to be imported, or is in ROOT distribution, and check if it can be persisted

  // Check if we already have the class (i.e. it is in the classToFile map)
  if (!doReplace && _c2fmap.find(tc->GetName())!=_c2fmap.end()) {
    oocxcoutD(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo(" << _wspace->GetName() << ") code of class " << tc->GetName() << " already imported, skipping" << endl ;
    return kTRUE ;
  }

  // Check if class is listed in a ROOTMAP file - if so we can skip it because it is in the root distribtion
  const char* mapEntry = gInterpreter->GetClassSharedLibs(tc->GetName()) ;
  if (mapEntry && strlen(mapEntry)>0) {
    oocxcoutD(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo(" << _wspace->GetName() << ") code of class " << tc->GetName() << " is in ROOT distribution, skipping " << endl ;
    return kTRUE ;
  }

  // Retrieve file names through ROOT TClass interface
  string implfile = tc->GetImplFileName() ;
  string declfile = tc->GetDeclFileName() ;

  // Check that file names are not empty
  if (implfile.empty() || declfile.empty()) {
    oocoutE(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo(" << _wspace->GetName() << ") ERROR: cannot retrieve code file names for class "
				   << tc->GetName() << " through ROOT TClass interface, unable to import code" << endl ;
    return kFALSE ;
  }

  // Check if header filename is found in ROOT distribution, if so, do not import class
  TString rootsys = gSystem->Getenv("ROOTSYS") ;
  if (TString(implfile.c_str()).Index(rootsys)>=0) {
    oocxcoutD(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo(" << _wspace->GetName() << ") code of class " << tc->GetName() << " is in ROOT distribution, skipping " << endl ;
    return kTRUE ;
  }

  // Require that class meets technical criteria to be persistable (i.e it has a default ctor)
  // (We also need a default ctor of abstract classes, but cannot check that through is interface
  //  as TClass::HasDefaultCtor only returns true for callable default ctors)
  if (!(tc->Property() & kIsAbstract) && !tc->HasDefaultConstructor()) {
    oocoutW(_wspace,ObjectHandling) << "RooWorkspace::autoImportClass(" << _wspace->GetName() << ") WARNING cannot import class "
				    << tc->GetName() << " : it cannot be persisted because it doesn't have a default constructor. Please fix " << endl ;
    return kFALSE ;
  }


  // *** PHASE 2 *** Check if declaration and implementation files can be located

  char *declpath = nullptr, *implpath = nullptr;

  // Check if header file can be found in specified location
  // If not, scan through list of 'class declaration' paths in RooWorkspace
  if (gSystem->AccessPathName(declfile.c_str())) {

    // Check list of additional declaration paths
    list<string>::iterator diter = RooWorkspace::_classDeclDirList.begin() ;

    while(diter!= RooWorkspace::_classDeclDirList.end()) {

      declpath = gSystem->ConcatFileName(diter->c_str(),declfile.c_str()) ;
      if (!gSystem->AccessPathName(declpath)) {
	// found declaration file
	break ;
      }
      // cleanup and continue ;
      delete [] declpath;
      declpath = nullptr;

      ++diter ;
    }

    // Header file cannot be found anywhere, warn user and abort operation
    if (!declpath) {
      oocoutW(_wspace,ObjectHandling) << "RooWorkspace::autoImportClass(" << _wspace->GetName() << ") WARNING Cannot access code of class "
				      << tc->GetName() << " because header file " << declfile << " is not found in current directory nor in $ROOTSYS" ;
      if (_classDeclDirList.size()>0) {
	ooccoutW(_wspace,ObjectHandling) << ", nor in the search path " ;
	diter = RooWorkspace::_classDeclDirList.begin() ;

	while(diter!= RooWorkspace::_classDeclDirList.end()) {

	  if (diter!=RooWorkspace::_classDeclDirList.begin()) {
	    ooccoutW(_wspace,ObjectHandling) << "," ;
	  }
	  ooccoutW(_wspace,ObjectHandling) << diter->c_str() ;
	  ++diter ;
	}
      }
      ooccoutW(_wspace,ObjectHandling) << ". To fix this problem, add the required directory to the search "
				       << "path using RooWorkspace::addClassDeclImportDir(const char* dir)" << endl ;

      return kFALSE ;
    }
  }


  // Check if implementation file can be found in specified location
  // If not, scan through list of 'class implementation' paths in RooWorkspace
  if (gSystem->AccessPathName(implfile.c_str())) {

    // Check list of additional declaration paths
    list<string>::iterator iiter = RooWorkspace::_classImplDirList.begin() ;

    while(iiter!= RooWorkspace::_classImplDirList.end()) {

      implpath = gSystem->ConcatFileName(iiter->c_str(),implfile.c_str()) ;
      if (!gSystem->AccessPathName(implpath)) {
	// found implementation file
	break ;
      }
      // cleanup and continue ;
      delete [] implpath;
      implpath = nullptr;

      ++iiter ;
    }

    // Implementation file cannot be found anywhere, warn user and abort operation
    if (!implpath) {
      oocoutW(_wspace,ObjectHandling) << "RooWorkspace::autoImportClass(" << _wspace->GetName() << ") WARNING Cannot access code of class "
				      << tc->GetName() << " because implementation file " << implfile << " is not found in current directory nor in $ROOTSYS" ;
      if (_classDeclDirList.size()>0) {
	ooccoutW(_wspace,ObjectHandling) << ", nor in the search path " ;
	iiter = RooWorkspace::_classImplDirList.begin() ;

	while(iiter!= RooWorkspace::_classImplDirList.end()) {

	  if (iiter!=RooWorkspace::_classImplDirList.begin()) {
	    ooccoutW(_wspace,ObjectHandling) << "," ;
	  }
	  ooccoutW(_wspace,ObjectHandling) << iiter->c_str() ;
	  ++iiter ;
	}
      }
      ooccoutW(_wspace,ObjectHandling) << ". To fix this problem add the required directory to the search "
				       << "path using RooWorkspace::addClassImplImportDir(const char* dir)" << endl;
      delete [] declpath;
      return kFALSE;
    }
  }

  char buf[64000];

  // *** Phase 3 *** Prepare to import code from files into STL string buffer
  //
  // Code storage is organized in two linked maps
  //
  // _fmap contains stl strings with code, indexed on declaration file name
  //
  // _c2fmap contains list of declaration file names and list of base classes
  //                  and is indexed on class name
  //
  // Phase 3 is skipped if fmap already contains an entry with given filebasename

  string declfilename = declpath?gSystem->BaseName(declpath):gSystem->BaseName(declfile.c_str()) ;

  // Split in base and extension
  int dotpos2 = strrchr(declfilename.c_str(),'.') - declfilename.c_str() ;
  string declfilebase = declfilename.substr(0,dotpos2) ;
  string declfileext = declfilename.substr(dotpos2+1) ;

  list<string> extraHeaders ;

  // If file has not beed stored yet, enter stl strings with implementation and declaration in file map
  if (_fmap.find(declfilebase) == _fmap.end()) {

    // Open declaration file
    fstream fdecl(declpath?declpath:declfile.c_str()) ;

    // Abort import if declaration file cannot be opened
    if (!fdecl) {
      oocoutE(_wspace,ObjectHandling) << "RooWorkspace::autoImportClass(" << _wspace->GetName()
				      << ") ERROR opening declaration file " <<  declfile << endl ;
      delete[] implpath;
      delete[] declpath;
      return kFALSE ;
    }

    oocoutI(_wspace,ObjectHandling) << "RooWorkspace::autoImportClass(" << _wspace->GetName()
				    << ") importing code of class " << tc->GetName()
				    << " from " << (implpath?implpath:implfile.c_str())
				    << " and " << (declpath?declpath:declfile.c_str()) << endl ;


    // Read entire file into an stl string
    string decl ;
    while(fdecl.getline(buf,1023)) {

      // Look for include state of self
      Bool_t processedInclude = kFALSE ;
      char* extincfile = 0 ;

      // Look for include of declaration file corresponding to this implementation file
      if (strstr(buf,"#include")) {
	// Process #include statements here
   char tmp[64000];
   strlcpy(tmp, buf, 64000);
   Bool_t stdinclude = strchr(buf, '<');
   strtok(tmp, " <\"");
   char *incfile = strtok(0, " <>\"");

   if (!stdinclude) {
      // check if it lives in $ROOTSYS/include
      TString hpath = gSystem->Getenv("ROOTSYS");
      hpath += "/include/";
      hpath += incfile;
      if (gSystem->AccessPathName(hpath.Data())) {
         oocoutI(_wspace, ObjectHandling) << "RooWorkspace::autoImportClass(" << _wspace->GetName()
                                          << ") scheduling include file " << incfile << " for import" << endl;
         extraHeaders.push_back(incfile);
         extincfile = incfile;
         processedInclude = kTRUE;
      }
   }
      }

      if (processedInclude) {
	decl += "// external include file below retrieved from workspace code storage\n" ;
	decl += Form("#include \"%s\"\n",extincfile) ;
      } else {
	decl += buf ;
	decl += '\n' ;
      }
    }

    // Open implementation file
    fstream fimpl(implpath?implpath:implfile.c_str()) ;

    // Abort import if implementation file cannot be opened
    if (!fimpl) {
      oocoutE(_wspace,ObjectHandling) << "RooWorkspace::autoImportClass(" << _wspace->GetName()
				      << ") ERROR opening implementation file " <<  implfile << endl ;
      delete[] implpath;
      delete[] declpath;
      return kFALSE ;
    }


    // Import entire implentation file into stl string
    string impl ;
    while(fimpl.getline(buf,1023)) {
      // Process #include statements here

      // Look for include state of self
      Bool_t foundSelfInclude=kFALSE ;
      Bool_t processedInclude = kFALSE ;
      char* extincfile = 0 ;

      // Look for include of declaration file corresponding to this implementation file
      if (strstr(buf,"#include")) {
	// Process #include statements here
   char tmp[64000];
   strlcpy(tmp, buf, 64000);
   Bool_t stdinclude = strchr(buf, '<');
   strtok(tmp, " <\"");
   char *incfile = strtok(0, " <>\"");

   if (strstr(incfile, declfilename.c_str())) {
      foundSelfInclude = kTRUE;
   }

   if (!stdinclude && !foundSelfInclude) {
      // check if it lives in $ROOTSYS/include
      TString hpath = gSystem->Getenv("ROOTSYS");
      hpath += "/include/";
      hpath += incfile;

      if (gSystem->AccessPathName(hpath.Data())) {
         oocoutI(_wspace, ObjectHandling) << "RooWorkspace::autoImportClass(" << _wspace->GetName()
                                          << ") scheduling include file " << incfile << " for import" << endl;
         extraHeaders.push_back(incfile);
         extincfile = incfile;
         processedInclude = kTRUE;
      }
   }
      }

      // Explicitly rewrite include of own declaration file to string
      // any directory prefixes, copy all other lines verbatim in stl string
      if (foundSelfInclude) {
	// If include of self is found, substitute original include
	// which may have directory structure with a plain include
	impl += "// class declaration include file below retrieved from workspace code storage\n" ;
	impl += Form("#include \"%s.%s\"\n",declfilebase.c_str(),declfileext.c_str()) ;
      } else if (processedInclude) {
	impl += "// external include file below retrieved from workspace code storage\n" ;
	impl += Form("#include \"%s\"\n",extincfile) ;
      } else {
	impl += buf ;
	impl += '\n' ;
      }
    }

    // Create entry in file map
    _fmap[declfilebase]._hfile = decl ;
    _fmap[declfilebase]._cxxfile = impl ;
    _fmap[declfilebase]._hext = declfileext ;

    // Process extra includes now
    for (list<string>::iterator ehiter = extraHeaders.begin() ; ehiter != extraHeaders.end() ; ++ehiter ) {
      if (_ehmap.find(*ehiter) == _ehmap.end()) {

	ExtraHeader eh ;
	eh._hname = ehiter->c_str() ;
	fstream fehdr(ehiter->c_str()) ;
	string ehimpl ;
	char buf2[1024] ;
	while(fehdr.getline(buf2,1023)) {

	  // Look for include of declaration file corresponding to this implementation file
	  if (strstr(buf2,"#include")) {
	    // Process #include statements here
       char tmp[64000];
       strlcpy(tmp, buf2, 64000);
       Bool_t stdinclude = strchr(buf, '<');
       strtok(tmp, " <\"");
       char *incfile = strtok(0, " <>\"");

       if (!stdinclude) {
          // check if it lives in $ROOTSYS/include
          TString hpath = gSystem->Getenv("ROOTSYS");
          hpath += "/include/";
          hpath += incfile;
          if (gSystem->AccessPathName(hpath.Data())) {
             oocoutI(_wspace, ObjectHandling) << "RooWorkspace::autoImportClass(" << _wspace->GetName()
                                              << ") scheduling recursive include file " << incfile << " for import"
                                              << endl;
             extraHeaders.push_back(incfile);
          }
       }
     }

     ehimpl += buf2;
     ehimpl += '\n';
   }
   eh._hfile = ehimpl.c_str();

   _ehmap[ehiter->c_str()] = eh;
      }
    }

  } else {

    // Inform that existing file entry is being recycled because it already contained class code
    oocoutI(_wspace,ObjectHandling) << "RooWorkspace::autoImportClass(" << _wspace->GetName()
				    << ") code of class " << tc->GetName()
				    << " was already imported from " << (implpath?implpath:implfile.c_str())
				    << " and " << (declpath?declpath:declfile.c_str()) << endl ;

  }


  // *** PHASE 4 *** Import stl strings with code into workspace
  //
  // If multiple classes are declared in a single code unit, there will be
  // multiple _c2fmap entries all pointing to the same _fmap entry.

  // Make list of all immediate base classes of this class
  TString baseNameList ;
  TList* bl = tc->GetListOfBases() ;
  TIterator* iter = bl->MakeIterator() ;
  TBaseClass* base ;
  list<TClass*> bases ;
  while((base=(TBaseClass*)iter->Next())) {
    if (baseNameList.Length()>0) {
      baseNameList += "," ;
    }
    baseNameList += base->GetClassPointer()->GetName() ;
    bases.push_back(base->GetClassPointer()) ;
  }

  // Map class name to above _fmap entries, along with list of base classes
  // in _c2fmap
  _c2fmap[tc->GetName()]._baseName = baseNameList ;
  _c2fmap[tc->GetName()]._fileBase = declfilebase ;

  // Recursive store all base classes.
  list<TClass*>::iterator biter = bases.begin() ;
  while(biter!=bases.end()) {
    autoImportClass(*biter,doReplace) ;
    ++biter ;
  }

  // Cleanup
  delete[] implpath;
  delete[] declpath;

  return kTRUE ;
}


////////////////////////////////////////////////////////////////////////////////
/// Create transient TDirectory representation of this workspace. This directory
/// will appear as a subdirectory of the directory that contains the workspace
/// and will have the name of the workspace suffixed with "Dir". The TDirectory
/// interface is read-only. Any attempt to insert objects into the workspace
/// directory representation will result in an error message. Note that some
/// ROOT object like TH1 automatically insert themselves into the current directory
/// when constructed. This will give error messages when done in a workspace
/// directory.

Bool_t RooWorkspace::makeDir()
{
  if (_dir) return kTRUE ;

  TString title= Form("TDirectory representation of RooWorkspace %s",GetName()) ;
  _dir = new WSDir(GetName(),title.Data(),this) ;

  TIterator* iter = componentIterator() ;
  RooAbsArg* darg ;
  while((darg=(RooAbsArg*)iter->Next())) {
    if (darg->IsA() != RooConstVar::Class()) {
      _dir->InternalAppend(darg) ;
    }
  }

  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Import a clone of a generic TObject into workspace generic object container. Imported
/// object can be retrieved by name through the obj() method. The object is cloned upon
/// importation and the input argument does not need to live beyond the import call
///
/// Returns kTRUE if an error has occurred.

Bool_t RooWorkspace::import(TObject& object, Bool_t replaceExisting)
{
  // First check if object with given name already exists
  TObject* oldObj = _genObjects.FindObject(object.GetName()) ;
  if (oldObj && !replaceExisting) {
    coutE(InputArguments) << "RooWorkspace::import(" << GetName() << ") generic object with name "
			  << object.GetName() << " is already in workspace and replaceExisting flag is set to false" << endl ;
    return kTRUE ;
  }

  // Grab the current state of the directory Auto-Add
  ROOT::DirAutoAdd_t func = object.IsA()->GetDirectoryAutoAdd();
  object.IsA()->SetDirectoryAutoAdd(0);
  Bool_t tmp = RooPlot::setAddDirectoryStatus(kFALSE) ;

  if (oldObj) {
    _genObjects.Replace(oldObj,object.Clone()) ;
    delete oldObj ;
  } else {
    _genObjects.Add(object.Clone()) ;
  }

  // Reset the state of the directory Auto-Add
  object.IsA()->SetDirectoryAutoAdd(func);
  RooPlot::setAddDirectoryStatus(tmp) ;

  return kFALSE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Import a clone of a generic TObject into workspace generic object container.
/// The imported object will be stored under the given alias name rather than its
/// own name. Imported object can be retrieved its alias name through the obj() method.
/// The object is cloned upon importation and the input argument does not need to live beyond the import call
/// This method is mostly useful for importing objects that do not have a settable name such as TMatrix
///
/// Returns kTRUE if an error has occurred.

Bool_t RooWorkspace::import(TObject& object, const char* aliasName, Bool_t replaceExisting)
{
  // First check if object with given name already exists
  TObject* oldObj = _genObjects.FindObject(object.GetName()) ;
  if (oldObj && !replaceExisting) {
    coutE(InputArguments) << "RooWorkspace::import(" << GetName() << ") generic object with name "
			  << object.GetName() << " is already in workspace and replaceExisting flag is set to false" << endl ;
    return kTRUE ;
  }

  TH1::AddDirectory(kFALSE) ;
  RooTObjWrap* wrapper = new RooTObjWrap(object.Clone()) ;
  TH1::AddDirectory(kTRUE) ;
  wrapper->setOwning(kTRUE) ;
  wrapper->SetName(aliasName) ;
  wrapper->SetTitle(aliasName) ;

  if (oldObj) {
    _genObjects.Replace(oldObj,wrapper) ;
    delete oldObj ;
  } else {
    _genObjects.Add(wrapper) ;
  }
  return kFALSE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Insert RooStudyManager module

Bool_t RooWorkspace::addStudy(RooAbsStudy& study)
{
  RooAbsStudy* clone = (RooAbsStudy*) study.Clone() ;
  _studyMods.Add(clone) ;
  return kFALSE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Remove all RooStudyManager modules

void RooWorkspace::clearStudies()
{
  _studyMods.Delete() ;
}




////////////////////////////////////////////////////////////////////////////////
/// Return any type of object (RooAbsArg, RooAbsData or generic object) with given name)

TObject* RooWorkspace::obj(const char* name) const
{
  // Try RooAbsArg first
  TObject* ret = arg(name) ;
  if (ret) return ret ;

  // Then try RooAbsData
  ret = data(name) ;
  if (ret) return ret ;

  // Finally try generic object store
  return genobj(name) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return generic object with given name

TObject* RooWorkspace::genobj(const char* name)  const
{
  // Find object by name
  TObject* gobj = _genObjects.FindObject(name) ;

  // Exit here if not found
  if (!gobj) return 0 ;

  // If found object is wrapper, return payload
  if (gobj->IsA()==RooTObjWrap::Class()) return ((RooTObjWrap*)gobj)->obj() ;

  return gobj ;
}



////////////////////////////////////////////////////////////////////////////////

Bool_t RooWorkspace::cd(const char* path)
{
  makeDir() ;
  return _dir->cd(path) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Save this current workspace into given file

Bool_t RooWorkspace::writeToFile(const char* fileName, Bool_t recreate)
{
  TFile f(fileName,recreate?"RECREATE":"UPDATE") ;
  Write() ;
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return instance to factory tool

RooFactoryWSTool& RooWorkspace::factory()
{
  if (_factory) {
    return *_factory;
  }
  cxcoutD(ObjectHandling) << "INFO: Creating RooFactoryWSTool associated with this workspace" << endl ;
  _factory = make_unique<RooFactoryWSTool>(*this);
  return *_factory;
}




////////////////////////////////////////////////////////////////////////////////
/// Short-hand function for `factory()->process(expr);`
///
/// \copydoc RooFactoryWSTool::process(const char*)
RooAbsArg* RooWorkspace::factory(const char* expr)
{
  return factory().process(expr) ;
}




////////////////////////////////////////////////////////////////////////////////
/// Print contents of the workspace

void RooWorkspace::Print(Option_t* opts) const
{
  Bool_t treeMode(kFALSE) ;
  Bool_t verbose(kFALSE);
  if (TString(opts).Contains("t")) {
    treeMode=kTRUE ;
  }
  if (TString(opts).Contains("v")) {
     verbose = kTRUE;
  }

  cout << endl << "RooWorkspace(" << GetName() << ") " << GetTitle() << " contents" << endl << endl  ;

  RooAbsArg* parg ;

  RooArgSet pdfSet ;
  RooArgSet funcSet ;
  RooArgSet varSet ;
  RooArgSet catfuncSet ;
  RooArgSet convResoSet ;
  RooArgSet resoSet ;


  // Split list of components in pdfs, functions and variables
  TIterator* iter = _allOwnedNodes.createIterator() ;
  while((parg=(RooAbsArg*)iter->Next())) {

    //---------------

    if (treeMode) {

      // In tree mode, only add nodes with no clients to the print lists

      if (parg->IsA()->InheritsFrom(RooAbsPdf::Class())) {
	if (!parg->hasClients()) {
	  pdfSet.add(*parg) ;
	}
      }

      if (parg->IsA()->InheritsFrom(RooAbsReal::Class()) &&
	  !parg->IsA()->InheritsFrom(RooAbsPdf::Class()) &&
	  !parg->IsA()->InheritsFrom(RooConstVar::Class()) &&
	  !parg->IsA()->InheritsFrom(RooRealVar::Class())) {
	if (!parg->hasClients()) {
	  funcSet.add(*parg) ;
	}
      }


      if (parg->IsA()->InheritsFrom(RooAbsCategory::Class()) &&
	  !parg->IsA()->InheritsFrom(RooCategory::Class())) {
	if (!parg->hasClients()) {
	  catfuncSet.add(*parg) ;
	}
      }

    } else {

      if (parg->IsA()->InheritsFrom(RooResolutionModel::Class())) {
	if (((RooResolutionModel*)parg)->isConvolved()) {
	  convResoSet.add(*parg) ;
	} else {
	  resoSet.add(*parg) ;
	}
      }

      if (parg->IsA()->InheritsFrom(RooAbsPdf::Class()) &&
	  !parg->IsA()->InheritsFrom(RooResolutionModel::Class())) {
	pdfSet.add(*parg) ;
      }

      if (parg->IsA()->InheritsFrom(RooAbsReal::Class()) &&
	  !parg->IsA()->InheritsFrom(RooAbsPdf::Class()) &&
	  !parg->IsA()->InheritsFrom(RooConstVar::Class()) &&
	  !parg->IsA()->InheritsFrom(RooRealVar::Class())) {
	funcSet.add(*parg) ;
      }

      if (parg->IsA()->InheritsFrom(RooAbsCategory::Class()) &&
	  !parg->IsA()->InheritsFrom(RooCategory::Class())) {
	catfuncSet.add(*parg) ;
      }
    }

    if (parg->IsA()->InheritsFrom(RooRealVar::Class())) {
      varSet.add(*parg) ;
    }

    if (parg->IsA()->InheritsFrom(RooCategory::Class())) {
      varSet.add(*parg) ;
    }

  }
  delete iter ;


  RooFit::MsgLevel oldLevel = RooMsgService::instance().globalKillBelow() ;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  if (varSet.getSize()>0) {
    varSet.sort() ;
    cout << "variables" << endl ;
    cout << "---------" << endl ;
    cout << varSet << endl ;
    cout << endl ;
  }

  if (pdfSet.getSize()>0) {
    cout << "p.d.f.s" << endl ;
    cout << "-------" << endl ;
    pdfSet.sort() ;
    iter = pdfSet.createIterator() ;
    while((parg=(RooAbsArg*)iter->Next())) {
      if (treeMode) {
	parg->printComponentTree() ;
      } else {
	parg->Print() ;
      }
    }
    delete iter ;
    cout << endl ;
  }

  if (!treeMode) {
    if (resoSet.getSize()>0) {
      cout << "analytical resolution models" << endl ;
      cout << "----------------------------" << endl ;
      resoSet.sort() ;
      iter = resoSet.createIterator() ;
      while((parg=(RooAbsArg*)iter->Next())) {
	parg->Print() ;
      }
      delete iter ;
      //     iter = convResoSet.createIterator() ;
      //     while((parg=(RooAbsArg*)iter->Next())) {
      //       parg->Print() ;
      //     }
      //     delete iter ;
      cout << endl ;
    }
  }

  if (funcSet.getSize()>0) {
    cout << "functions" << endl ;
    cout << "--------" << endl ;
    funcSet.sort() ;
    iter = funcSet.createIterator() ;
    while((parg=(RooAbsArg*)iter->Next())) {
      if (treeMode) {
	parg->printComponentTree() ;
      } else {
	parg->Print() ;
      }
    }
    delete iter ;
    cout << endl ;
  }

  if (catfuncSet.getSize()>0) {
    cout << "category functions" << endl ;
    cout << "------------------" << endl ;
    catfuncSet.sort() ;
    iter = catfuncSet.createIterator() ;
    while((parg=(RooAbsArg*)iter->Next())) {
      if (treeMode) {
	parg->printComponentTree() ;
      } else {
	parg->Print() ;
      }
    }
    delete iter ;
    cout << endl ;
  }

  if (_dataList.GetSize()>0) {
    cout << "datasets" << endl ;
    cout << "--------" << endl ;
    iter = _dataList.MakeIterator() ;
    RooAbsData* data2 ;
    while((data2=(RooAbsData*)iter->Next())) {
      cout << data2->IsA()->GetName() << "::" << data2->GetName() << *data2->get() << endl ;
    }
    delete iter ;
    cout << endl ;
  }

  if (_embeddedDataList.GetSize()>0) {
    cout << "embedded datasets (in pdfs and functions)" << endl ;
    cout << "-----------------------------------------" << endl ;
    iter = _embeddedDataList.MakeIterator() ;
    RooAbsData* data2 ;
    while((data2=(RooAbsData*)iter->Next())) {
      cout << data2->IsA()->GetName() << "::" << data2->GetName() << *data2->get() << endl ;
    }
    delete iter ;
    cout << endl ;
  }

  if (_snapshots.GetSize()>0) {
    cout << "parameter snapshots" << endl ;
    cout << "-------------------" << endl ;
    iter = _snapshots.MakeIterator() ;
    RooArgSet* snap ;
    while((snap=(RooArgSet*)iter->Next())) {
      cout << snap->GetName() << " = (" ;
      TIterator* aiter = snap->createIterator() ;
      RooAbsArg* a ;
      Bool_t first(kTRUE) ;
      while((a=(RooAbsArg*)aiter->Next())) {
	if (first) { first=kFALSE ; } else { cout << "," ; }
	cout << a->GetName() << "=" ;
	a->printValue(cout) ;
	if (a->isConstant()) {
	  cout << "[C]" ;
	}
      }
      cout << ")" << endl ;
      delete aiter ;
    }
    delete iter ;
    cout << endl ;
  }


  if (_namedSets.size()>0) {
    cout << "named sets" << endl ;
    cout << "----------" << endl ;
    for (map<string,RooArgSet>::const_iterator it = _namedSets.begin() ; it != _namedSets.end() ; ++it) {
       if (verbose || !TString(it->first.c_str()).BeginsWith("CACHE_")) {
          cout << it->first << ":" << it->second << endl;
       }
    }

    cout << endl ;
  }


  if (_genObjects.GetSize()>0) {
    cout << "generic objects" << endl ;
    cout << "---------------" << endl ;
    iter = _genObjects.MakeIterator() ;
    TObject* gobj ;
    while((gobj=(TObject*)iter->Next())) {
      if (gobj->IsA()==RooTObjWrap::Class()) {
	cout << ((RooTObjWrap*)gobj)->obj()->IsA()->GetName() << "::" << gobj->GetName() << endl ;
      } else {
	cout << gobj->IsA()->GetName() << "::" << gobj->GetName() << endl ;
      }
    }
    delete iter ;
    cout << endl ;

  }

  if (_studyMods.GetSize()>0) {
    cout << "study modules" << endl ;
    cout << "-------------" << endl ;
    iter = _studyMods.MakeIterator() ;
    TObject* smobj ;
    while((smobj=(TObject*)iter->Next())) {
      cout << smobj->IsA()->GetName() << "::" << smobj->GetName() << endl ;
    }
    delete iter ;
    cout << endl ;

  }

  if (_classes.listOfClassNames().size()>0) {
    cout << "embedded class code" << endl ;
    cout << "-------------------" << endl ;
    cout << _classes.listOfClassNames() << endl ;
    cout << endl ;
  }

  if (_eocache.size()>0) {
    cout << "embedded precalculated expensive components" << endl ;
    cout << "-------------------------------------------" << endl ;
    _eocache.print() ;
  }

  RooMsgService::instance().setGlobalKillBelow(oldLevel) ;

  return ;
}


////////////////////////////////////////////////////////////////////////////////
/// Custom streamer for the workspace. Stream contents of workspace
/// and code repository. When reading, read code repository first
/// and compile missing classes before proceeding with streaming
/// of workspace contents to avoid errors.

void RooWorkspace::CodeRepo::Streamer(TBuffer &R__b)
{
  typedef ::RooWorkspace::CodeRepo thisClass;

   // Stream an object of class RooWorkspace::CodeRepo.
   if (R__b.IsReading()) {

     UInt_t R__s, R__c;
     Version_t R__v =  R__b.ReadVersion(&R__s, &R__c);

     // Stream contents of ClassFiles map
     Int_t count(0) ;
     R__b >> count ;
     while(count--) {
       TString name ;
       name.Streamer(R__b) ;
       _fmap[name]._hext.Streamer(R__b) ;
       _fmap[name]._hfile.Streamer(R__b) ;
       _fmap[name]._cxxfile.Streamer(R__b) ;
     }

     // Stream contents of ClassRelInfo map
     count=0 ;
     R__b >> count ;
     while(count--) {
       TString name ;
       name.Streamer(R__b) ;
       _c2fmap[name]._baseName.Streamer(R__b) ;
       _c2fmap[name]._fileBase.Streamer(R__b) ;
     }

     if (R__v==2) {

       count=0 ;
       R__b >> count ;
       while(count--) {
	 TString name ;
	 name.Streamer(R__b) ;
	 _ehmap[name]._hname.Streamer(R__b) ;
	 _ehmap[name]._hfile.Streamer(R__b) ;
       }
     }

     R__b.CheckByteCount(R__s, R__c, thisClass::IsA());

     // Instantiate any classes that are not defined in current session
     _compiledOK = !compileClasses() ;

   } else {

     UInt_t R__c;
     R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);

     // Stream contents of ClassFiles map
     UInt_t count = _fmap.size() ;
     R__b << count ;
     map<TString,ClassFiles>::iterator iter = _fmap.begin() ;
     while(iter!=_fmap.end()) {
       TString key_copy(iter->first) ;
       key_copy.Streamer(R__b) ;
       iter->second._hext.Streamer(R__b) ;
       iter->second._hfile.Streamer(R__b);
       iter->second._cxxfile.Streamer(R__b);

       ++iter ;
     }

     // Stream contents of ClassRelInfo map
     count = _c2fmap.size() ;
     R__b << count ;
     map<TString,ClassRelInfo>::iterator iter2 = _c2fmap.begin() ;
     while(iter2!=_c2fmap.end()) {
       TString key_copy(iter2->first) ;
       key_copy.Streamer(R__b) ;
       iter2->second._baseName.Streamer(R__b) ;
       iter2->second._fileBase.Streamer(R__b);
       ++iter2 ;
     }

     // Stream contents of ExtraHeader map
     count = _ehmap.size() ;
     R__b << count ;
     map<TString,ExtraHeader>::iterator iter3 = _ehmap.begin() ;
     while(iter3!=_ehmap.end()) {
       TString key_copy(iter3->first) ;
       key_copy.Streamer(R__b) ;
       iter3->second._hname.Streamer(R__b) ;
       iter3->second._hfile.Streamer(R__b);
       ++iter3 ;
     }

     R__b.SetByteCount(R__c, kTRUE);

   }
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooWorkspace. This is a standard ROOT streamer for the
/// I/O part. This custom function exists to detach all external client links
/// from the payload prior to writing the payload so that these client links
/// are not persisted. (Client links occur if external function objects use
/// objects contained in the workspace as input)
/// After the actual writing, these client links are restored.

void RooWorkspace::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {

      R__b.ReadClassBuffer(RooWorkspace::Class(),this);

      // Perform any pass-2 schema evolution here
      RooFIter fiter = _allOwnedNodes.fwdIterator() ;
      RooAbsArg* node ;
      while((node=fiter.next())) {
	node->ioStreamerPass2() ;
      }
      RooAbsArg::ioStreamerPass2Finalize() ;

      // Make expensive object cache of all objects point to intermal copy.
      // Somehow this doesn't work OK automatically
      TIterator* iter = _allOwnedNodes.createIterator() ;
      while((node=(RooAbsArg*)iter->Next())) {
	node->setExpensiveObjectCache(_eocache) ;
   node->setWorkspace(*this);
   if (node->IsA()->InheritsFrom(RooAbsOptTestStatistic::Class())) {
      RooAbsOptTestStatistic *tmp = (RooAbsOptTestStatistic *)node;
      if (tmp->isSealed() && tmp->sealNotice() && strlen(tmp->sealNotice()) > 0) {
         cout << "RooWorkspace::Streamer(" << GetName() << ") " << node->IsA()->GetName() << "::" << node->GetName()
              << " : " << tmp->sealNotice() << endl;
      }
   }
      }
      delete iter ;


   } else {

     // Make lists of external clients of WS objects, and remove those links temporarily

     map<RooAbsArg*,vector<RooAbsArg *> > extClients, extValueClients, extShapeClients ;

     TIterator* iter = _allOwnedNodes.createIterator() ;
     RooAbsArg* tmparg ;
     while((tmparg=(RooAbsArg*)iter->Next())) {

       // Loop over client list of this arg
       std::vector<RooAbsArg *> clientsTmp{tmparg->_clientList.begin(), tmparg->_clientList.end()};
       for (auto client : clientsTmp) {
         if (!_allOwnedNodes.containsInstance(*client)) {

           const auto refCount = tmparg->_clientList.refCount(client);
           auto& bufferVec = extClients[tmparg];

           bufferVec.insert(bufferVec.end(), refCount, client);
           tmparg->_clientList.Remove(client, true);
         }
       }

       // Loop over value client list of this arg
       clientsTmp.assign(tmparg->_clientListValue.begin(), tmparg->_clientListValue.end());
       for (auto vclient : clientsTmp) {
         if (!_allOwnedNodes.containsInstance(*vclient)) {
           cxcoutD(ObjectHandling) << "RooWorkspace::Streamer(" << GetName() << ") element " << tmparg->GetName()
				       << " has external value client link to " << vclient << " (" << vclient->GetName() << ") with ref count " << tmparg->_clientListValue.refCount(vclient) << endl ;

           const auto refCount = tmparg->_clientListValue.refCount(vclient);
           auto& bufferVec = extValueClients[tmparg];

           bufferVec.insert(bufferVec.end(), refCount, vclient);
           tmparg->_clientListValue.Remove(vclient, true);
         }
       }

       // Loop over shape client list of this arg
       clientsTmp.assign(tmparg->_clientListShape.begin(), tmparg->_clientListShape.end());
       for (auto sclient : clientsTmp) {
         if (!_allOwnedNodes.containsInstance(*sclient)) {
           cxcoutD(ObjectHandling) << "RooWorkspace::Streamer(" << GetName() << ") element " << tmparg->GetName()
				         << " has external shape client link to " << sclient << " (" << sclient->GetName() << ") with ref count " << tmparg->_clientListShape.refCount(sclient) << endl ;

           const auto refCount = tmparg->_clientListShape.refCount(sclient);
           auto& bufferVec = extShapeClients[tmparg];

           bufferVec.insert(bufferVec.end(), refCount, sclient);
           tmparg->_clientListShape.Remove(sclient, true);
         }
       }

     }
     delete iter ;

     R__b.WriteClassBuffer(RooWorkspace::Class(),this);

     // Reinstate clients here


     for (auto iterx : extClients) {
       for (auto client : iterx.second) {
         iterx.first->_clientList.Add(client);
       }
     }

     for (auto iterx : extValueClients) {
       for (auto client : iterx.second) {
         iterx.first->_clientListValue.Add(client);
       }
     }

     for (auto iterx : extShapeClients) {
       for (auto client : iterx.second) {
         iterx.first->_clientListShape.Add(client);
       }
     }

   }
}




////////////////////////////////////////////////////////////////////////////////
/// Return STL string with last of class names contained in the code repository

std::string RooWorkspace::CodeRepo::listOfClassNames() const
{
  string ret ;
  map<TString,ClassRelInfo>::const_iterator iter = _c2fmap.begin() ;
  while(iter!=_c2fmap.end()) {
    if (ret.size()>0) {
      ret += ", " ;
    }
    ret += iter->first ;
    ++iter ;
  }

  return ret ;
}

namespace {
UInt_t crc32(const char* data, ULong_t sz, UInt_t crc)
{
  // update CRC32 with new data

  // use precomputed table, rather than computing it on the fly
  static const UInt_t crctab[256] = { 0x00000000,
    0x04c11db7, 0x09823b6e, 0x0d4326d9, 0x130476dc, 0x17c56b6b,
    0x1a864db2, 0x1e475005, 0x2608edb8, 0x22c9f00f, 0x2f8ad6d6,
    0x2b4bcb61, 0x350c9b64, 0x31cd86d3, 0x3c8ea00a, 0x384fbdbd,
    0x4c11db70, 0x48d0c6c7, 0x4593e01e, 0x4152fda9, 0x5f15adac,
    0x5bd4b01b, 0x569796c2, 0x52568b75, 0x6a1936c8, 0x6ed82b7f,
    0x639b0da6, 0x675a1011, 0x791d4014, 0x7ddc5da3, 0x709f7b7a,
    0x745e66cd, 0x9823b6e0, 0x9ce2ab57, 0x91a18d8e, 0x95609039,
    0x8b27c03c, 0x8fe6dd8b, 0x82a5fb52, 0x8664e6e5, 0xbe2b5b58,
    0xbaea46ef, 0xb7a96036, 0xb3687d81, 0xad2f2d84, 0xa9ee3033,
    0xa4ad16ea, 0xa06c0b5d, 0xd4326d90, 0xd0f37027, 0xddb056fe,
    0xd9714b49, 0xc7361b4c, 0xc3f706fb, 0xceb42022, 0xca753d95,
    0xf23a8028, 0xf6fb9d9f, 0xfbb8bb46, 0xff79a6f1, 0xe13ef6f4,
    0xe5ffeb43, 0xe8bccd9a, 0xec7dd02d, 0x34867077, 0x30476dc0,
    0x3d044b19, 0x39c556ae, 0x278206ab, 0x23431b1c, 0x2e003dc5,
    0x2ac12072, 0x128e9dcf, 0x164f8078, 0x1b0ca6a1, 0x1fcdbb16,
    0x018aeb13, 0x054bf6a4, 0x0808d07d, 0x0cc9cdca, 0x7897ab07,
    0x7c56b6b0, 0x71159069, 0x75d48dde, 0x6b93dddb, 0x6f52c06c,
    0x6211e6b5, 0x66d0fb02, 0x5e9f46bf, 0x5a5e5b08, 0x571d7dd1,
    0x53dc6066, 0x4d9b3063, 0x495a2dd4, 0x44190b0d, 0x40d816ba,
    0xaca5c697, 0xa864db20, 0xa527fdf9, 0xa1e6e04e, 0xbfa1b04b,
    0xbb60adfc, 0xb6238b25, 0xb2e29692, 0x8aad2b2f, 0x8e6c3698,
    0x832f1041, 0x87ee0df6, 0x99a95df3, 0x9d684044, 0x902b669d,
    0x94ea7b2a, 0xe0b41de7, 0xe4750050, 0xe9362689, 0xedf73b3e,
    0xf3b06b3b, 0xf771768c, 0xfa325055, 0xfef34de2, 0xc6bcf05f,
    0xc27dede8, 0xcf3ecb31, 0xcbffd686, 0xd5b88683, 0xd1799b34,
    0xdc3abded, 0xd8fba05a, 0x690ce0ee, 0x6dcdfd59, 0x608edb80,
    0x644fc637, 0x7a089632, 0x7ec98b85, 0x738aad5c, 0x774bb0eb,
    0x4f040d56, 0x4bc510e1, 0x46863638, 0x42472b8f, 0x5c007b8a,
    0x58c1663d, 0x558240e4, 0x51435d53, 0x251d3b9e, 0x21dc2629,
    0x2c9f00f0, 0x285e1d47, 0x36194d42, 0x32d850f5, 0x3f9b762c,
    0x3b5a6b9b, 0x0315d626, 0x07d4cb91, 0x0a97ed48, 0x0e56f0ff,
    0x1011a0fa, 0x14d0bd4d, 0x19939b94, 0x1d528623, 0xf12f560e,
    0xf5ee4bb9, 0xf8ad6d60, 0xfc6c70d7, 0xe22b20d2, 0xe6ea3d65,
    0xeba91bbc, 0xef68060b, 0xd727bbb6, 0xd3e6a601, 0xdea580d8,
    0xda649d6f, 0xc423cd6a, 0xc0e2d0dd, 0xcda1f604, 0xc960ebb3,
    0xbd3e8d7e, 0xb9ff90c9, 0xb4bcb610, 0xb07daba7, 0xae3afba2,
    0xaafbe615, 0xa7b8c0cc, 0xa379dd7b, 0x9b3660c6, 0x9ff77d71,
    0x92b45ba8, 0x9675461f, 0x8832161a, 0x8cf30bad, 0x81b02d74,
    0x857130c3, 0x5d8a9099, 0x594b8d2e, 0x5408abf7, 0x50c9b640,
    0x4e8ee645, 0x4a4ffbf2, 0x470cdd2b, 0x43cdc09c, 0x7b827d21,
    0x7f436096, 0x7200464f, 0x76c15bf8, 0x68860bfd, 0x6c47164a,
    0x61043093, 0x65c52d24, 0x119b4be9, 0x155a565e, 0x18197087,
    0x1cd86d30, 0x029f3d35, 0x065e2082, 0x0b1d065b, 0x0fdc1bec,
    0x3793a651, 0x3352bbe6, 0x3e119d3f, 0x3ad08088, 0x2497d08d,
    0x2056cd3a, 0x2d15ebe3, 0x29d4f654, 0xc5a92679, 0xc1683bce,
    0xcc2b1d17, 0xc8ea00a0, 0xd6ad50a5, 0xd26c4d12, 0xdf2f6bcb,
    0xdbee767c, 0xe3a1cbc1, 0xe760d676, 0xea23f0af, 0xeee2ed18,
    0xf0a5bd1d, 0xf464a0aa, 0xf9278673, 0xfde69bc4, 0x89b8fd09,
    0x8d79e0be, 0x803ac667, 0x84fbdbd0, 0x9abc8bd5, 0x9e7d9662,
    0x933eb0bb, 0x97ffad0c, 0xafb010b1, 0xab710d06, 0xa6322bdf,
    0xa2f33668, 0xbcb4666d, 0xb8757bda, 0xb5365d03, 0xb1f740b4
  };

  crc = ~crc;
  while (sz--) crc = (crc << 8) ^ UInt_t(*data++) ^ crctab[crc >> 24];

  return ~crc;
}

UInt_t crc32(const char* data)
{
  // Calculate crc32 checksum on given string
  unsigned long sz = strlen(data);
  switch (strlen(data)) {
    case 0:
      return 0;
    case 1:
      return data[0];
    case 2:
      return (data[0] << 8) | data[1];
    case 3:
      return (data[0] << 16) | (data[1] << 8) | data[2];
    case 4:
      return (data[0] << 24) | (data[1] << 16) | (data[2] << 8) | data[3];
    default:
      return crc32(data + 4, sz - 4, (data[0] << 24) | (data[1] << 16) |
      (data[2] << 8) | data[3]);
  }
}

}

////////////////////////////////////////////////////////////////////////////////
/// For all classes in the workspace for which no class definition is
/// found in the ROOT class table extract source code stored in code
/// repository into temporary directory set by
/// setClassFileExportDir(), compile classes and link them with
/// current ROOT session. If a compilation error occurs print
/// instructions for user how to fix errors and recover workspace and
/// abort import procedure.

Bool_t RooWorkspace::CodeRepo::compileClasses()
{
  Bool_t haveDir=kFALSE ;

  // Retrieve name of directory in which to export code files
  string dirName = Form(_classFileExportDir.c_str(),_wspace->uuid().AsString(),_wspace->GetName()) ;

  Bool_t writeExtraHeaders(kFALSE) ;

  // Process all class entries in repository
  map<TString,ClassRelInfo>::iterator iter = _c2fmap.begin() ;
  while(iter!=_c2fmap.end()) {

    oocxcoutD(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() now processing class " << iter->first.Data() << endl ;

    // If class is already known, don't load
    if (gClassTable->GetDict(iter->first.Data())) {
      oocoutI(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() Embedded class "
				      << iter->first << " already in ROOT class table, skipping" << endl ;
      ++iter ;
      continue ;
    }

    // Check that export directory exists
    if (!haveDir) {

      // If not, make local directory to extract files
      if (!gSystem->AccessPathName(dirName.c_str())) {
	oocoutI(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() reusing code export directory " << dirName.c_str()
					<< " to extract coded embedded in workspace" << endl ;
      } else {
	if (gSystem->MakeDirectory(dirName.c_str())==0) {
	  oocoutI(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() creating code export directory " << dirName.c_str()
					  << " to extract coded embedded in workspace" << endl ;
	} else {
	  oocoutE(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() ERROR creating code export directory " << dirName.c_str()
					  << " to extract coded embedded in workspace" << endl ;
	  return kFALSE ;
	}
      }
      haveDir=kTRUE ;

    }

    // First write any extra header files
    if (!writeExtraHeaders) {
      writeExtraHeaders = kTRUE ;

      map<TString,ExtraHeader>::iterator eiter = _ehmap.begin() ;
      while(eiter!=_ehmap.end()) {

	// Check if identical declaration file (header) is already written
	Bool_t needEHWrite=kTRUE ;
	string fdname = Form("%s/%s",dirName.c_str(),eiter->second._hname.Data()) ;
	ifstream ifdecl(fdname.c_str()) ;
	if (ifdecl) {
	  TString contents ;
     char buf[64000];
     while (ifdecl.getline(buf, 64000)) {
        contents += buf;
        contents += "\n";
     }
     UInt_t crcFile = crc32(contents.Data());
     UInt_t crcWS = crc32(eiter->second._hfile.Data());
     needEHWrite = (crcFile != crcWS);
   }

   // Write declaration file if required
   if (needEHWrite) {
      oocoutI(_wspace, ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() Extracting extra header file "
                                       << fdname << endl;

      // Extra headers may contain non-existing path - create first to be sure
      gSystem->MakeDirectory(gSystem->GetDirName(fdname.c_str()));

      ofstream fdecl(fdname.c_str());
      if (!fdecl) {
         oocoutE(_wspace, ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() ERROR opening file " << fdname
                                          << " for writing" << endl;
         return kFALSE;
      }
      fdecl << eiter->second._hfile.Data();
      fdecl.close();
   }
   ++eiter;
      }
    }


    // Navigate from class to file
    ClassFiles& cfinfo = _fmap[iter->second._fileBase] ;

    oocxcoutD(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() now processing file with base " << iter->second._fileBase << endl ;

    // If file is already processed, skip to next class
    if (cfinfo._extracted) {
      oocxcoutD(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() file with base name " << iter->second._fileBase
					 << " has already been extracted, skipping to next class" << endl ;
      continue ;
    }

    // Check if identical declaration file (header) is already written
    Bool_t needDeclWrite=kTRUE ;
    string fdname = Form("%s/%s.%s",dirName.c_str(),iter->second._fileBase.Data(),cfinfo._hext.Data()) ;
    ifstream ifdecl(fdname.c_str()) ;
    if (ifdecl) {
      TString contents ;
      char buf[64000];
      while (ifdecl.getline(buf, 64000)) {
         contents += buf;
         contents += "\n";
      }
      UInt_t crcFile = crc32(contents.Data()) ;
      UInt_t crcWS   = crc32(cfinfo._hfile.Data()) ;
      needDeclWrite = (crcFile!=crcWS) ;
    }

    // Write declaration file if required
    if (needDeclWrite) {
      oocoutI(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() Extracting declaration code of class " << iter->first << ", file " << fdname << endl ;
      ofstream fdecl(fdname.c_str()) ;
      if (!fdecl) {
	oocoutE(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() ERROR opening file "
					<< fdname << " for writing" << endl ;
	return kFALSE ;
      }
      fdecl << cfinfo._hfile ;
      fdecl.close() ;
    }

    // Check if identical implementation file is already written
    Bool_t needImplWrite=kTRUE ;
    string finame = Form("%s/%s.cxx",dirName.c_str(),iter->second._fileBase.Data()) ;
    ifstream ifimpl(finame.c_str()) ;
    if (ifimpl) {
      TString contents ;
      char buf[64000];
      while (ifimpl.getline(buf, 64000)) {
         contents += buf;
         contents += "\n";
      }
      UInt_t crcFile = crc32(contents.Data()) ;
      UInt_t crcWS   = crc32(cfinfo._cxxfile.Data()) ;
      needImplWrite = (crcFile!=crcWS) ;
    }

    // Write implementation file if required
    if (needImplWrite) {
      oocoutI(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() Extracting implementation code of class " << iter->first << ", file " << finame << endl ;
      ofstream fimpl(finame.c_str()) ;
      if (!fimpl) {
	oocoutE(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() ERROR opening file"
					<< finame << " for writing" << endl ;
	return kFALSE ;
      }
      fimpl << cfinfo._cxxfile ;
      fimpl.close() ;
    }

    // Mark this file as extracted
    cfinfo._extracted = kTRUE ;
    oocxcoutD(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() marking code unit  " << iter->second._fileBase << " as extracted" << endl ;

    // Compile class
    oocoutI(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() Compiling code unit " << iter->second._fileBase.Data() << " to define class " << iter->first << endl ;
    Bool_t ok = gSystem->CompileMacro(finame.c_str(),"k") ;

    if (!ok) {
      oocoutE(_wspace,ObjectHandling) << "RooWorkspace::CodeRepo::compileClasses() ERROR compiling class " << iter->first.Data() << ", to fix this you can do the following: " << endl
				      << "  1) Fix extracted source code files in directory " << dirName.c_str() << "/" << endl
				      << "  2) In clean ROOT session compiled fixed classes by hand using '.x " << dirName.c_str() << "/ClassName.cxx+'" << endl
				      << "  3) Reopen file with RooWorkspace with broken source code in UPDATE mode. Access RooWorkspace to force loading of class" << endl
				      << "     Broken instances in workspace will _not_ be compiled, instead precompiled fixed instances will be used." << endl
				      << "  4) Reimport fixed code in workspace using 'RooWorkspace::importClassCode(\"*\",kTRUE)' method, Write() updated workspace to file and close file" << endl
				      << "  5) Reopen file in clean ROOT session to confirm that problems are fixed" << endl ;
	return kFALSE ;
    }

    ++iter ;
  }

  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Internal access to TDirectory append method

void RooWorkspace::WSDir::InternalAppend(TObject* obj)
{
  TDirectory::Append(obj,kFALSE) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Overload TDirectory interface method to prohibit insertion of objects in read-only directory workspace representation

void RooWorkspace::WSDir::Add(TObject* obj,Bool_t)
{
  if (dynamic_cast<RooAbsArg*>(obj) || dynamic_cast<RooAbsData*>(obj)) {
    coutE(ObjectHandling) << "RooWorkspace::WSDir::Add(" << GetName() << ") ERROR: Directory is read-only representation of a RooWorkspace, use RooWorkspace::import() to add objects" << endl ;
  } else {
    InternalAppend(obj) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Overload TDirectory interface method to prohibit insertion of objects in read-only directory workspace representation

void RooWorkspace::WSDir::Append(TObject* obj,Bool_t)
{
  if (dynamic_cast<RooAbsArg*>(obj) || dynamic_cast<RooAbsData*>(obj)) {
    coutE(ObjectHandling) << "RooWorkspace::WSDir::Add(" << GetName() << ") ERROR: Directory is read-only representation of a RooWorkspace, use RooWorkspace::import() to add objects" << endl ;
  } else {
    InternalAppend(obj) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Activate export of workspace symbols to CINT in a namespace with given name. If no name
/// is given the namespace will have the same name as the workspace

void RooWorkspace::exportToCint(const char* nsname)
{
  // If export is already active, do nothing
  if (_doExport) {
    coutE(ObjectHandling) << "RooWorkspace::exportToCint(" << GetName() << ") WARNING: repeated calls to exportToCint() have no effect" << endl ;
    return ;
  }

  // Set flag so that future import to workspace are automatically exported to CINT
  _doExport = kTRUE ;

  // If no name is provided choose name of workspace
  if (!nsname) nsname = GetName() ;
  _exportNSName = nsname ;

  coutI(ObjectHandling) << "RooWorkspace::exportToCint(" << GetName()
			<< ") INFO: references to all objects in this workspace will be created in CINT in 'namespace " << _exportNSName << "'" << endl ;

  // Export present contents of workspace to CINT
  TIterator* iter = _allOwnedNodes.createIterator() ;
  TObject* wobj ;
  while((wobj=iter->Next())) {
    exportObj(wobj) ;
  }
  delete iter ;
  iter = _dataList.MakeIterator() ;
  while((wobj=iter->Next())) {
    exportObj(wobj) ;
  }
  delete iter ;
}


////////////////////////////////////////////////////////////////////////////////
/// Export reference to given workspace object to CINT

void RooWorkspace::exportObj(TObject* wobj)
{
  // Do nothing if export flag is not set
  if (!_doExport) return ;

  // Do not export RooConstVars
  if (wobj->IsA() == RooConstVar::Class()) {
    return ;
  }


  // Determine if object name is a valid C++ identifier name

  // Do not export objects that have names that are not valid C++ identifiers
  if (!isValidCPPID(wobj->GetName())) {
    cxcoutD(ObjectHandling) << "RooWorkspace::exportObj(" << GetName() << ") INFO: Workspace object name " << wobj->GetName() << " is not a valid C++ identifier and is not exported to CINT" << endl ;
    return ;
  }

  // Declare correctly typed reference to object in CINT in the namespace associated with this workspace
  string cintExpr = Form("namespace %s { %s& %s = *(%s *)0x%lx ; }",_exportNSName.c_str(),wobj->IsA()->GetName(),wobj->GetName(),wobj->IsA()->GetName(),(ULong_t)wobj) ;
  gROOT->ProcessLine(cintExpr.c_str()) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return true if given name is a valid C++ identifier name

Bool_t RooWorkspace::isValidCPPID(const char* name)
{
  string oname(name) ;
  if (isdigit(oname[0])) {
    return kFALSE ;
  } else {
    for (UInt_t i=0 ; i<oname.size() ; i++) {
      char c = oname[i] ;
      if (!isalnum(c) && (c!='_')) {
	return kFALSE ;
      }
    }
  }
  return kTRUE ;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete exported reference in CINT namespace

void RooWorkspace::unExport()
{
   char buf[64000];
   TIterator *iter = _allOwnedNodes.createIterator();
   TObject *wobj;
   while ((wobj = iter->Next())) {
      if (isValidCPPID(wobj->GetName())) {
         strlcpy(buf, Form("%s::%s", _exportNSName.c_str(), wobj->GetName()), 64000);
         gInterpreter->DeleteVariable(buf);
      }
  }
  delete iter ;
}

////////////////////////////////////////////////////////////////////////////////
/// If one of the TObject we have a referenced to is deleted, remove the
/// reference.

void RooWorkspace::RecursiveRemove(TObject *removedObj)
{
   _dataList.RecursiveRemove(removedObj);
   if (removedObj == _dir) _dir = nullptr;

   _allOwnedNodes.RecursiveRemove(removedObj); // RooArgSet

   _dataList.RecursiveRemove(removedObj);
   _embeddedDataList.RecursiveRemove(removedObj);
   _views.RecursiveRemove(removedObj);
   _snapshots.RecursiveRemove(removedObj);
   _genObjects.RecursiveRemove(removedObj);
   _studyMods.RecursiveRemove(removedObj);

   for(auto c : _namedSets) {
      c.second.RecursiveRemove(removedObj);
   }

   _eocache.RecursiveRemove(removedObj); // RooExpensiveObjectCache
}
