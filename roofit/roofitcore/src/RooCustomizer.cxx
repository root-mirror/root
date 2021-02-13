/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
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

/**
 * \class RooCustomizer
 *
 * RooCustomizer is a factory class to produce clones
 * of a prototype composite PDF object with the same structure but
 * different leaf servers (parameters or dependents).
 *
 * RooCustomizer supports two kinds of modifications:
 *
 * - replaceArg(leaf_arg, repl_arg):
 * Replaces each occurence of leaf_arg with repl_arg in the composite pdf.
 *
 * - splitArg(split_arg):
 * Build multiple clones of the same prototype. Each
 * occurrence of split_arg is replaced with a clone of split_arg
 * named split_arg_[MCstate], where [MCstate] is the name of the
 * 'master category state' that indexes the clones to be built.
 *
 *
 * ### Example: Change the decay constant of an exponential for each run
 *
 * Splitting is particularly useful when building simultaneous fits to
 * subsets of the data sample with different background properties.
 * In such a case, the user builds a single prototype PDF representing
 * the structure of the signal and background and splits the dataset
 * into categories with different background properties. Using
 * RooCustomizer a PDF for each subfit can be constructed from the
 * prototype that has same structure and signal parameters, but
 * different instances of the background parameters: e.g.
 * ```
 *     ...
 *     RooExponential bg("bg","background",x,alpha) ;
 *     RooGaussian sig("sig","signal",x,mean,sigma) ;
 *     RooAddPdf pdf("pdf","pdf",sig,bg,sigfrac) ;
 *
 *     RooDataSet data("data","dataset",RooArgSet(x,runblock),...)
 *
 *     RooCategory runblock("runblock","run block") ;
 *     runblock.defineType("run1") ;
 *     runblock.defineType("run2") ;
 *
 *     RooArgSet splitLeafs;
 *     RooCustomizer cust(pdf,runblock,splitLeafs);
 *     cust.splitArg(alpha,runblock);
 *
 *     RooAbsPdf* pdf_run1 = cust.build("run1") ;
 *     RooAbsPdf* pdf_run2 = cust.build("run2") ;
 *
 *     RooSimultaneous simpdf("simpdf","simpdf",RooArgSet(*pdf_run1,*pdf_run2))
 * ```
 * If the master category state is a super category, leafs may be split
 * by any subset of that master category. E.g. if the master category
 * is 'A x B', leafs may be split by A, B or AxB.
 *
 * In addition to replacing leaf nodes, RooCustomizer clones all branch
 * nodes that depend directly or indirectly on modified leaf nodes, so
 * that the input pdf is untouched by each build operation.
 *
 * The customizer owns all the branch nodes including the returned top
 * level node, so the customizer should live as longs as the cloned
 * composites are needed.
 *
 * Any leaf nodes that are created by the customizer will be put into
 * the leaf list that is passed into the customizers constructor (splitLeafs in
 * the above example. The list owner is responsible for deleting these leaf
 * nodes after the customizer is deleted.
 *
 *
 * ## Advanced techniques
 *
 * ### Reuse nodes to customise a different PDF
 * By default, the customizer clones the prototype leaf node when splitting a leaf,
 * but the user can feed pre-defined split leafs in leaf list. These leafs
 * must have the name <split_leaf>_<splitcat_label> to be picked up. The list
 * of pre-supplied leafs may be partial, any missing split leafs will be auto
 * generated.
 *
 * Another common construction is to have two prototype PDFs, each to be customized
 * by a separate customizer instance, that share parameters. To ensure that
 * the customized clones also share their respective split leafs, i.e.
 * ```
 *   PDF1(x,y, A) and PDF2(z, A) ---> PDF1_run1(x,y, A_run1) and PDF2_run1(x,y, A_run1)
 *                                    PDF1_run2(x,y, A_run2) and PDF2_run2(x,y, A_run2)
 * ```
 * feed the same split leaf list into both customizers. In that case, the second customizer
 * will pick up the split leafs instantiated by the first customizer and the link between
 * the two PDFs is retained.
 *
 * ### Customising with pre-defined leafs
 * If leaf nodes are provided in the sets, the customiser will use them. This is a complete
 * example that customises the `yield` parameter, and splits (automatically clones) the
 * mean of the Gaussian. This is a short version of the tutorial rf514_RooCustomizer.C.
 * ```
 *  RooRealVar E("Energy","Energy",0,3000);
 *
 *  RooRealVar meanG("meanG","meanG", peak[1]);
 *  RooRealVar fwhm("fwhm", "fwhm", 5/(2*Sqrt(2*Log(2))));
 *  RooGaussian gauss("gauss", "gauss", E, meanG, fwhm);
 *
 *  RooPolynomial linear("linear","linear",E,RooArgList());
 *
 *  RooRealVar yieldSig("yieldSig", "yieldSig", 1, 0, 1.E4);
 *  RooRealVar yieldBkg("yieldBkg", "yieldBkg", 1, 0, 1.E4);
 *
 *  RooAddPdf model("model","model",
 *      RooArgList(gauss,linear),
 *      RooArgList(yieldSig, yieldBkg));
 *
 *  RooCategory sample("sample","sample");
 *  sample.defineType("BBG1m2T");
 *  sample.defineType("BBG2m2T");
 *
 *
 *  RooArgSet customisedLeafs;
 *  RooArgSet allLeafs;
 *
 *  RooRealVar mass("M", "M", 1, 0, 12000);
 *  RooFormulaVar yield1("yieldSig_BBG1m2T","sigy1","M/3.360779",mass);
 *  RooFormulaVar yield2("yieldSig_BBG2m2T","sigy2","M/2",mass);
 *  allLeafs.add(yield1);
 *  allLeafs.add(yield2);
 *
 *
 *  RooCustomizer cust(model, sample, customisedLeafs, &allLeafs);
 *  cust.splitArg(yieldSig, sample);
 *  cust.splitArg(meanG, sample);
 *
 *  auto pdf1 = cust.build("BBG1m2T");
 *  auto pdf2 = cust.build("BBG2m2T");
 * ```
*/


#include "RooCustomizer.h"

#include "RooAbsCategoryLValue.h"
#include "RooAbsCategory.h"
#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooMsgService.h"
#include "RooHelpers.h"

#include <iostream>
#include "strtok.h"
#include "strlcpy.h"

#include "RooWorkspace.h"
#include "RooGlobalFunc.h"
#include "RooConstVar.h"
#include "RooRealConstant.h"


#ifndef _WIN32
#include <strings.h>
#endif


using namespace std;

ClassImp(RooCustomizer); 
;

namespace {

static Int_t init();

Int_t dummy = init() ;

static Int_t init()
{
  RooFactoryWSTool::IFace* iface = new RooCustomizer::CustIFace ;
  RooFactoryWSTool::registerSpecial("EDIT",iface) ;
  (void) dummy;
  return 0 ;
}

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with a prototype and masterCat index category.
/// Customizers created by this constructor offer both the
/// replaceArg() and splitArg() functionality.
/// \param[in] pdf Proto PDF to be customised.
/// \param[in] masterCat Category to be used for splitting.
/// \param[in,out] splitLeafs All nodes created in
/// the customisation process are added to this set.
/// The user can provide nodes that are *taken*
/// from the set if they have a name that matches `<parameterNameToBeReplaced>_<category>`.
/// \note The set needs to own its contents if they are user-provided.
/// Use *e.g.*
/// ```
///  RooArgSet customisedLeafs;
///  auto yield1 = new RooFormulaVar("yieldSig_BBG1m2T","sigy1","M/3.360779",mass);
///  customisedLeafs.addOwned(*yield1);
/// ```
/// \param[in,out] splitLeafsAll All leafs that are used when customising are collected here.
/// If this set already contains leaves, they will be used for customising if the names match
/// as above.
/// 

RooCustomizer::RooCustomizer(const RooAbsArg& pdf, const RooAbsCategoryLValue& masterCat,
    RooArgSet& splitLeafs, RooArgSet* splitLeafsAll) :
  TNamed(pdf.GetName(),pdf.GetTitle()),
  _sterile(kFALSE),
  _owning(kTRUE),
  _masterPdf((RooAbsArg*)&pdf), 
  _masterCat((RooAbsCategoryLValue*)&masterCat), 
  _masterBranchList("masterBranchList"), 
  _masterLeafList("masterLeafList"), 
  _internalCloneBranchList("cloneBranchList"),
  _cloneNodeListAll(splitLeafsAll),
  _cloneNodeListOwned(&splitLeafs)
{
  _cloneBranchList = &_internalCloneBranchList ;

  initialize() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Sterile Constructor. Customizers created by this constructor
/// offer only the replace() method. The supplied 'name' is used as
/// suffix for any cloned branch nodes

RooCustomizer::RooCustomizer(const RooAbsArg& pdf, const char* name) :
  TNamed(pdf.GetName(),pdf.GetTitle()),
  _sterile(kTRUE), 
  _owning(kFALSE),
  _name(name),
  _masterPdf((RooAbsArg*)&pdf), 
  _masterCat(0), 
  _masterBranchList("masterBranchList"), 
  _masterLeafList("masterLeafList"), 
  _internalCloneBranchList("cloneBranchList"),
  _cloneNodeListAll(0),
  _cloneNodeListOwned(0)
{
  _cloneBranchList = &_internalCloneBranchList ;

  initialize() ;
}




////////////////////////////////////////////////////////////////////////////////
/// Initialize the customizer

void RooCustomizer::initialize() 
{
  _masterPdf->leafNodeServerList(&_masterLeafList) ;
  _masterPdf->branchNodeServerList(&_masterBranchList) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooCustomizer::~RooCustomizer() 
{

}


  

////////////////////////////////////////////////////////////////////////////////
/// Split all arguments in 'set' into individualized clones for each
/// defined state of 'splitCat'. The 'splitCats' category must be
/// subset of or equal to the master category supplied in the
/// customizer constructor.
///
/// Splitting is only available on customizers created with a master index category

void RooCustomizer::splitArgs(const RooArgSet& set, const RooAbsCategory& splitCat) 
{
  if (_sterile) {
    coutE(InputArguments) << "RooCustomizer::splitArgs(" << _name 
			  << ") ERROR cannot set spitting rules on this sterile customizer" << endl ;
    return ;
  }

  for (auto arg : set) {
    splitArg(*arg,splitCat) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Split all argument 'arg' into individualized clones for each
/// defined state of 'splitCat'. The 'splitCats' category must be
/// subset of or equal to the master category supplied in the
/// customizer constructor.
///
/// Splitting is only available on customizers created with a master index category

void RooCustomizer::splitArg(const RooAbsArg& arg, const RooAbsCategory& splitCat) 
{
  if (_splitArgList.FindObject(arg.GetName())) {
    coutE(InputArguments) << "RooCustomizer(" << GetName() << ") ERROR: multiple splitting rules defined for " 
			  << arg.GetName() << " only using first rule" << endl ;
    return ;
  }

  if (_sterile) {
    coutE(InputArguments) << "RooCustomizer::splitArg(" << _name 
	 << ") ERROR cannot set spitting rules on this sterile customizer" << endl ;
    return ;
  }

  _splitArgList.Add((RooAbsArg*)&arg) ;
  _splitCatList.Add((RooAbsCategory*)&splitCat) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Replace any occurence of arg 'orig' with arg 'subst'

void RooCustomizer::replaceArg(const RooAbsArg& orig, const RooAbsArg& subst) 
{
  if (_replaceArgList.FindObject(orig.GetName())) {
    coutE(InputArguments) << "RooCustomizer(" << GetName() << ") ERROR: multiple replacement rules defined for " 
	 << orig.GetName() << " only using first rule" << endl ;
    return ;
  }

  _replaceArgList.Add((RooAbsArg*)&orig) ;
  _replaceSubList.Add((RooAbsArg*)&subst) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Build a clone of the prototype executing all registered 'replace' rules.
/// If verbose is set, a message is printed for each leaf or branch node
/// modification. The returned head node owns all cloned branch nodes
/// that were created in the cloning process.

RooAbsArg* RooCustomizer::build(Bool_t verbose) 
{
  // Execute build
  RooAbsArg* ret =  doBuild(_name.Length()>0?_name.Data():0,verbose) ;

  // Make root object own all cloned nodes

  // First make list of all objects that were created
  RooArgSet allOwned ;
  if (_cloneNodeListOwned) {
    allOwned.add(*_cloneNodeListOwned) ;
  }
  allOwned.add(*_cloneBranchList) ;

  // Remove head node from list
  allOwned.remove(*ret) ;

  // If list with owned objects is not empty, assign
  // head node as owner
  if (allOwned.getSize()>0) {
    ret->addOwnedComponents(allOwned) ;
  }

  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Build a clone of the prototype executing all registered 'replace'
/// rules and 'split' rules for the masterCat state named
/// 'masterCatState'.  If verbose is set a message is printed for
/// each leaf or branch node modification. The returned composite arg
/// is owned by the customizer.  This function cannot be called on
/// customizer build with the sterile constructor.

RooAbsArg* RooCustomizer::build(const char* masterCatState, Bool_t verbose) 
{
  if (_sterile) {
    coutE(InputArguments) << "RooCustomizer::build(" << _name 
			  << ") ERROR cannot use leaf spitting build() on this sterile customizer" << endl ;
    return 0 ;
  }

  // Set masterCat to given state
  if (_masterCat->setLabel(masterCatState)) {
    coutE(InputArguments) << "RooCustomizer::build(" << _masterPdf->GetName() << "): ERROR label '" << masterCatState 
			  << "' not defined for master splitting category " << _masterCat->GetName() << endl ;
    return 0 ;
  }

  return doBuild(masterCatState,verbose) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Back-end implementation of the p.d.f building functionality

RooAbsArg* RooCustomizer::doBuild(const char* masterCatState, Bool_t verbose) 
{
  // Find nodes that must be split according to provided description, Clone nodes, change their names
  RooArgSet masterNodesToBeSplit("masterNodesToBeSplit") ;
  RooArgSet masterNodesToBeReplaced("masterNodesToBeReplaced") ;
  RooArgSet masterReplacementNodes("masterReplacementNodes") ;
  RooArgSet clonedMasterNodes("clonedMasterNodes") ;


  RooArgSet nodeList(_masterLeafList) ;
  nodeList.add(_masterBranchList) ;

  //   cout << "loop over " << nodeList.getSize() << " nodes" << endl ;
  for (auto node : nodeList) {
    RooAbsArg* theSplitArg = !_sterile?(RooAbsArg*) _splitArgList.FindObject(node->GetName()):0 ;
    if (theSplitArg) {
      RooAbsCategory* splitCat = (RooAbsCategory*) _splitCatList.At(_splitArgList.IndexOf(theSplitArg)) ;
      if (verbose) {
	coutI(ObjectHandling) << "RooCustomizer::build(" << _masterPdf->GetName() 
			      << "): tree node " << node->GetName() << " is split by category " << splitCat->GetName() << endl ;
      }
      
      TString newName(node->GetName()) ;
      if (masterCatState) {
	newName.Append("_") ;
	newName.Append(splitCat->getCurrentLabel()) ;	
      }

      // Check if this node instance already exists
      RooAbsArg* specNode = _cloneNodeListAll ? _cloneNodeListAll->find(newName) : _cloneNodeListOwned->find(newName) ;
      if (specNode) {

	// Copy instance to one-time use list for this build
	clonedMasterNodes.add(*specNode) ;
	if (verbose) {
	  coutI(ObjectHandling) << "RooCustomizer::build(" << _masterPdf->GetName() 
				<< ") Adding existing node specialization " << newName << " to clonedMasterNodes" << endl ;
	}

	// Affix attribute with old name to clone to support name changing server redirect
	TString nameAttrib("ORIGNAME:") ;
	nameAttrib.Append(node->GetName()) ;
	specNode->setAttribute(nameAttrib) ;
	
	if (!specNode->getStringAttribute("origName")) {
	  specNode->setStringAttribute("origName",node->GetName()) ;
	}



      } else {

	if (node->isDerived()) {
	  coutW(ObjectHandling) << "RooCustomizer::build(" << _masterPdf->GetName() 
				<< "): WARNING: branch node " << node->GetName() << " is split but has no pre-defined specializations" << endl ;
	}

	TString newTitle(node->GetTitle()) ;
	newTitle.Append(" (") ;
	newTitle.Append(splitCat->getCurrentLabel()) ;
	newTitle.Append(")") ;
      
	// Create a new clone
	RooAbsArg* clone = (RooAbsArg*) node->Clone(newName.Data()) ;
	clone->setStringAttribute("factory_tag",0) ;
	clone->SetTitle(newTitle) ;

	// Affix attribute with old name to clone to support name changing server redirect
	TString nameAttrib("ORIGNAME:") ;
	nameAttrib.Append(node->GetName()) ;
	clone->setAttribute(nameAttrib) ;

	if (!clone->getStringAttribute("origName")) {
	  clone->setStringAttribute("origName",node->GetName()) ;
	}

	// Add to one-time use list and life-time use list
	clonedMasterNodes.add(*clone) ;
	if (_owning) {
	  _cloneNodeListOwned->addOwned(*clone) ;	
	} else {
	  _cloneNodeListOwned->add(*clone) ;	
	}
	if (_cloneNodeListAll) {
	  _cloneNodeListAll->add(*clone) ;	
	}	
      }
      masterNodesToBeSplit.add(*node) ;     
    }

    RooAbsArg* ReplaceArg = (RooAbsArg*) _replaceArgList.FindObject(node->GetName()) ;
    if (ReplaceArg) {
      RooAbsArg* substArg = (RooAbsArg*) _replaceSubList.At(_replaceArgList.IndexOf(ReplaceArg)) ;
      if (verbose) {
	coutI(ObjectHandling) << "RooCustomizer::build(" << _masterPdf->GetName() 
			      << "): tree node " << node->GetName() << " will be replaced by " << substArg->GetName() << endl ;
      }

      // Affix attribute with old name to support name changing server redirect
      TString nameAttrib("ORIGNAME:") ;
      nameAttrib.Append(node->GetName()) ;
      substArg->setAttribute(nameAttrib) ;

      // Add to list
      masterNodesToBeReplaced.add(*node) ;
      masterReplacementNodes.add(*substArg) ;
    }
  }

  // Find branches that are affected by splitting and must be cloned
  RooArgSet masterBranchesToBeCloned("masterBranchesToBeCloned") ;
  for (auto branch : _masterBranchList) {
    
    // If branch is split itself, don't handle here
    if (masterNodesToBeSplit.find(branch->GetName())) {
      if (verbose) {
	coutI(ObjectHandling) << "RooCustomizer::build(" << _masterPdf->GetName() << ") Branch node " << branch->GetName() << " is already split" << endl ;
      }
      continue ;
    }
    if (masterNodesToBeReplaced.find(branch->GetName())) {
      if (verbose) {
	coutI(ObjectHandling) << "RooCustomizer::build(" << _masterPdf->GetName() << ") Branch node " << branch->GetName() << " is already replaced" << endl ;
      }
      continue ;
    }

    if (branch->dependsOn(masterNodesToBeSplit)) {
      if (verbose) {
	coutI(ObjectHandling) << "RooCustomizer::build(" << _masterPdf->GetName() << ") Branch node " 
			      << branch->IsA()->GetName() << "::" << branch->GetName() << " cloned: depends on a split parameter" << endl ;
      }
      masterBranchesToBeCloned.add(*branch) ;
    } else if (branch->dependsOn(masterNodesToBeReplaced)) {
      if (verbose) {
	coutI(ObjectHandling) << "RooCustomizer::build(" << _masterPdf->GetName() << ") Branch node " 
			      << branch->IsA()->GetName() << "::" << branch->GetName() << " cloned: depends on a replaced parameter" << endl ;
      }
      masterBranchesToBeCloned.add(*branch) ;
    }
  }

  // Clone branches, changes their names 
  RooAbsArg* cloneTopPdf = 0;
  RooArgSet clonedMasterBranches("clonedMasterBranches") ;

  for (auto branch : masterBranchesToBeCloned) {
    TString newName(branch->GetName()) ;
    if (masterCatState) {
      newName.Append("_") ;
      newName.Append(masterCatState) ;
    }

    // Affix attribute with old name to clone to support name changing server redirect
    RooAbsArg* clone = (RooAbsArg*) branch->Clone(newName.Data()) ;
    clone->setStringAttribute("factory_tag",0) ;
    TString nameAttrib("ORIGNAME:") ;
    nameAttrib.Append(branch->GetName()) ;
    clone->setAttribute(nameAttrib) ;

    if (!clone->getStringAttribute("origName")) {
      clone->setStringAttribute("origName",branch->GetName()) ;
    }
    
    clonedMasterBranches.add(*clone) ;      

    // Save pointer to clone of top-level pdf
    if (branch==_masterPdf) cloneTopPdf=(RooAbsArg*)clone ;
  }

  if (_owning) {
    _cloneBranchList->addOwned(clonedMasterBranches) ;
  } else {
    _cloneBranchList->add(clonedMasterBranches) ;
  }
    
  // Reconnect cloned branches to each other and to cloned nodess
  for (auto branch : clonedMasterBranches) {
    branch->redirectServers(clonedMasterBranches,kFALSE,kTRUE) ;
    branch->redirectServers(clonedMasterNodes,kFALSE,kTRUE) ;
    branch->redirectServers(masterReplacementNodes,kFALSE,kTRUE) ;
  }

  return cloneTopPdf ? cloneTopPdf : _masterPdf ;
}


////////////////////////////////////////////////////////////////////////////////
/// Print name of customizer

void RooCustomizer::printName(ostream& os) const 
{
  os << GetName() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Print title of customizer

void RooCustomizer::printTitle(ostream& os) const 
{
  os << GetTitle() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Print class name of customizer

void RooCustomizer::printClassName(ostream& os) const 
{
  os << IsA()->GetName() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Print arguments of customizer, i.e. input p.d.f and input master category (if any)

void RooCustomizer::printArgs(ostream& os) const 
{
  os << "[ masterPdf=" << _masterPdf->GetName() ;
  if (_masterCat) {
    os << " masterCat=" << _masterCat->GetName() ;
  }
  os << " ]" ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print customizer configuration details

void RooCustomizer::printMultiline(ostream& os, Int_t /*content*/, Bool_t /*verbose*/, TString indent) const
{
  os << indent << "RooCustomizer for " << _masterPdf->GetName() << (_sterile?" (sterile)":"") << endl ;

  Int_t i, nsplit = _splitArgList.GetSize() ;
  if (nsplit>0) {
    os << indent << "  Splitting rules:" << endl ;
    for (i=0 ; i<nsplit ; i++) {
      os << indent << "   " << _splitArgList.At(i)->GetName() << " is split by " << _splitCatList.At(i)->GetName() << endl ;
    }
  }

  Int_t nrepl = _replaceArgList.GetSize() ;
  if (nrepl>0) {
    os << indent << "  Replacement rules:" << endl ;
    for (i=0 ; i<nrepl ; i++) {
      os << indent << "   " << _replaceSubList.At(i)->GetName() << " replaces " << _replaceArgList.At(i)->GetName() << endl ;
    }
  }
  
  return ;
}



////////////////////////////////////////////////////////////////////////////////
/// Install the input RooArgSet as container in which all cloned branches
/// will be stored

void RooCustomizer::setCloneBranchSet(RooArgSet& cloneBranchSet) 
{
  _cloneBranchList = &cloneBranchSet ;
  _cloneBranchList->useHashMapForFind(true);
}




////////////////////////////////////////////////////////////////////////////////

std::string RooCustomizer::CustIFace::create(RooFactoryWSTool& ft, const char* typeName, const char* instanceName, std::vector<std::string> args) 
{
  // Check number of arguments
  if (args.size()<2) {
    throw string(Form("RooCustomizer::CustIFace::create() ERROR: expect at least 2 arguments for EDIT: the input object and at least one $Replace() rule")) ;
  }

  if (string(typeName)!="EDIT") {
    throw string(Form("RooCustomizer::CustIFace::create() ERROR: unknown type requested: %s",typeName)) ;
  }

  // Check that first arg exists as RooAbsArg
  RooAbsArg* arg = ft.ws().arg(args[0].c_str()) ;
  if (!arg) {
    throw string(Form("RooCustomizer::CustIFace::create() ERROR: input RooAbsArg %s does not exist",args[0].c_str())) ;
  }

  // If name of new object is same as original, execute in sterile mode (i.e no suffixes attached), and rename original nodes in workspace upon import
  if (args[0]==instanceName) {
    instanceName=0 ;
  }

  // Create a customizer
  RooCustomizer cust(*arg,instanceName) ;
  
  for (unsigned int i=1 ; i<args.size() ; i++) {
    char buf[1024] ;
    strlcpy(buf,args[i].c_str(),1024) ;
    char* sep = strchr(buf,'=') ;
    if (!sep) {
      throw string(Form("RooCustomizer::CustIFace::create() ERROR: unknown argument: %s, expect form orig=subst",args[i].c_str())) ;
    }
    *sep = 0 ;    
    RooAbsArg* orig = ft.ws().arg(buf) ;
    RooAbsArg* subst(0) ;
    if (string(sep+1).find("$REMOVE")==0) {

      // Create a removal dummy ;
      subst = &RooRealConstant::removalDummy() ;

      // If removal instructed was annotated with target node, encode these in removal dummy
      char* sep2 = strchr(sep+1,'(') ;
      if (sep2) {
	char buf2[1024] ;
	strlcpy(buf2,sep2+1,1024) ;
	char* saveptr ;
	char* tok = R__STRTOK_R(buf2,",)",&saveptr) ;
	while(tok) {
	  //cout << "$REMOVE is restricted to " << tok << endl ;
	  subst->setAttribute(Form("REMOVE_FROM_%s",tok)) ;
	  tok = R__STRTOK_R(0,",)",&saveptr) ;
	}
      } else {
	// Otherwise mark as universal removal node
	subst->setAttribute("REMOVE_ALL") ;
      }

    } else {
      subst = ft.ws().arg(sep+1) ;
    }
//     if (!orig) {
//       throw string(Form("RooCustomizer::CustIFace::create() ERROR: $Replace() input RooAbsArg %s does not exist",buf)) ;
//     }
//     if (!subst) {
//       throw string(Form("RooCustomizer::CustIFace::create() ERROR: $Replace() replacement RooAbsArg %s does not exist",sep+1)) ;
//     }
    if (orig && subst) {
      cust.replaceArg(*orig,*subst) ;
    } else {
      oocoutW((TObject*)0,ObjectHandling) << "RooCustomizer::CustIFace::create() WARNING: input or replacement of a replacement operation not found, operation ignored"<< endl ;
    }
  }

  // Build the desired edited object
  RooAbsArg* targ = cust.build(kFALSE)  ;
  if (!targ) {
    throw string(Form("RooCustomizer::CustIFace::create() ERROR in customizer build, object %snot created",instanceName)) ;
  }

  // Import the object into the workspace
  if (instanceName) {
    // Set the desired name of the top level node
    targ->SetName(instanceName) ;
    // Now import everything. What we didn't touch gets recycled, everything else was cloned here:
    ft.ws().import(cust.cloneBranchList(), RooFit::Silence(true), RooFit::RecycleConflictNodes(true),    RooFit::NoRecursion(false));
  } else {
    ft.ws().import(cust.cloneBranchList(), RooFit::Silence(true), RooFit::RenameConflictNodes("orig",1), RooFit::NoRecursion(true));
  }

  return string(instanceName?instanceName:targ->GetName()) ;
}
