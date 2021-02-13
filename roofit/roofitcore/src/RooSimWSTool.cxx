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

//////////////////////////////////////////////////////////////////////////////
/// \class RooSimWSTool
/// The RooSimWSTool is a tool operating on RooWorkspace objects that
/// can clone PDFs into a series of variations that are joined together
/// into a RooSimultanous PDF.
///
/// ## Splitting a single PDF
/// The simplest use case is to take a workspace PDF as prototype and
/// "split" a parameter of that PDF into two specialized parameters
/// depending on a category in the dataset.
///
/// For example, given a Gaussian
/// PDF \f$ G(x \,|\, m,s) \f$ we want to construct a \f$ G_a(x \,|\, m_a,s) \f$ and a \f$ G_b(x \,|\, m_b,s) \f$
/// with different mean parameters to be fit to a dataset with observables
/// \f$ (x,c) \f$ where \f$ c \f$ is a category with states 'a' and 'b'.
///
/// Using RooSimWSTool, one can create a simultaneous PDF from \f$ G_a \f$ and \f$ G_b \f$
/// from \f$ G \f$ with the following commands:
/// ```
///   RooSimWSTool wst(wspace);
///   wst.build("G_sim", "G", SplitParam("m","c"));
/// ```
///
/// #### Splitting using a product category
/// From this simple example one can go to builds of arbitrary complexity
/// by specifying multiple SplitParam arguments on multiple parameters
/// involving multiple splitting categories. Splits can also be performed
/// in the product of multiple categories, *i.e.*,
/// ```
///   wst.build("G_sim", "G", SplitParam("m","c,d"));
/// ```
/// splits the parameter \f$ m \f$ in the product of the states of \f$ c \f$ and
/// \f$ d \f$.
///
/// #### Constrained split
/// Another possibility
/// is the "constrained" split, which clones the parameter for all but one state
/// and inserts a formula specialization in a chosen state that evaluates
/// to \f$ 1 - \sum_i(a_i) \f$ where \f$ a_i \f$ are all other specializations. For example,
/// given a category \f$ c \f$ with the states `"A","B","C","D"`, the specification
/// ```
///     SplitParamConstrained("m","c","D")
/// ```
/// will create the parameters \f$ m_A,m_B,m_C \f$ and a formula expression \f$ m_D \f$
/// that evaluates to \f$ (1-(m_A+m_B+m_C)) \f$. Constrained splits can also be
/// specified in the product of categories. In that case, the name of the
/// remainder state follows the syntax `"{State1;State2}"`, where `State1` and
/// `State2` are the state names of the two spitting categories.
///
/// ## Splitting multiple PDFs
/// The examples so far deal with a single prototype PDF. It is also
/// possible to build with multiple prototype PDFs by specifying a
/// mapping between the prototype to use and the names of states of
/// a "master" splitting category. To specify these configurations,
/// an intermediate `MultiBuildConfig` must be composed with all
/// the necessary specifications. This, for example,
/// ```
///  RooSimWSTool::MultiBuildConfig mbc("mc");
///  mbc.addPdf("I","G",SplitParam("m,s","c"));
///  mbc.addPdf("II,III","F",SplitParam("a","c,d"));
/// ```
/// configures a build with two prototype PDFs \f$ G \f$ and \f$ F \f$.
/// Prototype \f$ G \f$ is used for state `"I"` of the master split category
/// `mc` and prototype \f$ F \f$ is used for states `"II"` and `"III"` of the
/// master split category `mc`. Furthermore, the parameters \f$ m,s \f$ of
/// prototype \f$ G \f$ are split
/// in category \f$ c \f$ while the parameter \f$ a \f$ of prototype \f$ F \f$ is split in
/// the product of the categories \f$ c \f$ and \f$ d \f$. The actual build is then
/// performed by passing the build configuration to RooSimWSTool, *e.g.*,
/// ```
///  wst.build("MASTER", mbc);
/// ```
/// By default, a specialisation is built for each permutation of
/// states of the splitting categories that are used. It is possible
/// to restrict the building of specialised PDFs to a subset of states
/// by adding a restriction on the number of states to build as follows:
/// ```
///  mbc.restrictBuild("c","A,B");
/// ```
/// The restrictBuild method can be called multiple times, but at most
/// once for each splitting category in use. For simple builds with a single
/// prototype, restriction can be specified with a Restrict() argument
/// on the build command line.
///


#include "RooFit.h"
#include "RooSimWSTool.h"

#include "RooFitLegacy/RooCatTypeLegacy.h"
#include "RooMsgService.h"
#include "RooCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooSuperCategory.h"
#include "RooCustomizer.h"
#include "RooMultiCategory.h"
#include "RooSimultaneous.h"
#include "RooGlobalFunc.h"
#include "RooFracRemainder.h"
#include "RooFactoryWSTool.h"

ClassImp(RooSimWSTool); 
ClassImp(RooSimWSTool::BuildConfig); 
ClassImp(RooSimWSTool::MultiBuildConfig); 
ClassImp(RooSimWSTool::SplitRule); 
ClassImp(RooSimWSTool::ObjBuildConfig); 
ClassImp(RooSimWSTool::ObjSplitRule); 
;

using namespace std ;

namespace {

static Int_t init();

Int_t dummy = init() ;

static Int_t init()
{
  RooFactoryWSTool::IFace* iface = new RooSimWSTool::SimWSIFace ;
  RooFactoryWSTool::registerSpecial("SIMCLONE",iface) ;
  RooFactoryWSTool::registerSpecial("MSIMCLONE",iface) ;
  (void) dummy;
  return 0 ;
}

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor of SimWSTool on given workspace. All input is taken from the workspace
/// All output is stored in the workspace

RooSimWSTool::RooSimWSTool(RooWorkspace& ws) : _ws(&ws) 
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooSimWSTool::~RooSimWSTool() 
{
}



////////////////////////////////////////////////////////////////////////////////
/// Build a RooSimultaneous PDF with name simPdfName from cloning specializations of protytpe PDF protoPdfName.
/// <table>
/// <tr><th> Optional Arguments <th> Effect
/// <tr><td> SplitParam(varname, catname)                   <td> Split parameter(s) with given name(s) in category(s) with given names
/// <tr><td> SplitParam(var, cat)                           <td> Split given parameter(s) in givem category(s)
/// <tr><td> SplitParamConstrained(vname, cname, remainder) <td> Make constrained split in parameter(s) with given name(s) in category(s) with given names
///                                                   putting remainder fraction formula in state with name "remainder"
/// <tr><td> SplitParamConstrained(var,cat,remainder)       <td> Make constrained split in parameter(s) with given name(s) in category(s) with given names
///                                                   putting remainder fraction formula in state with name "remainder"
/// <tr><td> Restrict(catName,stateNameList)                <td> Restrict build by only considered listed state names of category with given name

RooSimultaneous* RooSimWSTool::build(const char* simPdfName, const char* protoPdfName, const RooCmdArg& arg1,const RooCmdArg& arg2,
					const RooCmdArg& arg3,const RooCmdArg& arg4, const RooCmdArg& arg5,const RooCmdArg& arg6)
{
  BuildConfig bc(protoPdfName,arg1,arg2,arg3,arg4,arg5,arg6) ;
  return build(simPdfName,bc) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Build a RooSimultaneous PDF with name simPdfName from cloning specializations of protytpe PDF protoPdfName.
/// Use the provided BuildConfig or MultiBuildConfig object to configure the build

RooSimultaneous* RooSimWSTool::build(const char* simPdfName,BuildConfig& bc, Bool_t verbose) 
{
  ObjBuildConfig* obc = validateConfig(bc) ;
  if (!obc) return 0 ;
  
  if (verbose) {
    obc->print() ;
  }
  
  RooSimultaneous* ret =  executeBuild(simPdfName,*obc,verbose) ;

  delete obc ;
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Validate build configuration. If not syntax errors or missing objects are found,
/// return an ObjBuildConfig in which all names are replaced with object pointers.

RooSimWSTool::ObjBuildConfig* RooSimWSTool::validateConfig(BuildConfig& bc)
{
  // Create empty object version of build config
  ObjBuildConfig* obc = new ObjBuildConfig ;

  if (bc._masterCatName.length()>0) {
    obc->_masterCat = _ws->cat(bc._masterCatName.c_str()) ;
    if (!obc->_masterCat) {
      coutE(ObjectHandling) << "RooSimWSTool::build(" << GetName() << ") ERROR: associated workspace " << _ws->GetName() 
			    << " does not contain a category named " << bc._masterCatName 
			    << " that was designated as master index category in the build configuration" << endl ;
      delete obc ;
      return 0 ;
    }
  } else {
    obc->_masterCat = 0 ;
  }
  
  map<string,SplitRule>::iterator pdfiter ;
  // Check that we have the p.d.f.s
  for (pdfiter = bc._pdfmap.begin() ; pdfiter != bc._pdfmap.end() ; ++pdfiter) {    

    // Check that p.d.f exists
    RooAbsPdf* pdf = _ws->pdf(pdfiter->second.GetName()) ;
    if (!pdf) {
      coutE(ObjectHandling) << "RooSimWSTool::build(" << GetName() << ") ERROR: associated workspace " << _ws->GetName() 
			    << " does not contain a pdf named " << pdfiter->second.GetName() << endl ;
      delete obc ;
      return 0 ;
    }      

    // Create empty object version of split rule set
    ObjSplitRule osr ;

    // Convert names of parameters and splitting categories to objects in workspace, fill object split rule
    SplitRule& sr = pdfiter->second ;

    map<string, pair<list<string>,string> >::iterator pariter ;
    for (pariter=sr._paramSplitMap.begin() ; pariter!=sr._paramSplitMap.end() ; ++pariter) {
      
      // Check that variable with given name exists in workspace
      RooAbsArg* farg = _ws->fundArg(pariter->first.c_str()) ;
      if (!farg) {
	coutE(ObjectHandling) << "RooSimWSTool::build(" << GetName() << ") ERROR: associated workspace " << _ws->GetName() 
			      << " does not contain a variable named " << pariter->first.c_str() 
			      << " as specified in splitting rule of parameter " << pariter->first << " of p.d.f " << pdf << endl ;
	delete obc ;
	return 0 ;
      } 
      
      // Check that given variable is indeed related to given p.d.f
      if (!pdf->dependsOn(*farg)) {
	coutE(ObjectHandling) << "RooSimWSTool::build(" << GetName() << ") ERROR: specified parameter " << pariter->first 
			      << " in split is not function of p.d.f " << pdf->GetName() << endl ;
	delete obc ;
	return 0 ;
      }
      
      
      RooArgSet splitCatSet ;
      list<string>::iterator catiter ;
      for (catiter = pariter->second.first.begin() ; catiter!=pariter->second.first.end() ; ++catiter) {
	RooAbsCategory* cat = _ws->catfunc(catiter->c_str()) ;
	if (!cat) {
	  coutE(ObjectHandling) << "RooSimWSTool::build(" << GetName() << ") ERROR: associated workspace " << _ws->GetName() 
				<< " does not contain a category named " << catiter->c_str() 
				<< " as specified in splitting rule of parameter " << pariter->first << " of p.d.f " << pdf << endl ;
	  delete obc ;
	  return 0 ;
	}
	splitCatSet.add(*cat) ;
      }      

      // Check if composite splitCatSet does not contain category functions that depend on other categories used in the same split
      TIterator* iter = splitCatSet.createIterator() ;
      RooAbsArg* arg ;
      while((arg=(RooAbsArg*)iter->Next())) {
	RooArgSet tmp(splitCatSet) ;
	tmp.remove(*arg) ;
	if (arg->dependsOnValue(tmp)) {
	  coutE(InputArguments) << "RooSimWSTool::build(" << GetName() << ") ERROR: Ill defined split: splitting category function " << arg->GetName() 
				<< " used in composite split " << splitCatSet << " of parameter " << farg->GetName() << " of pdf " << pdf->GetName() 
				<< " depends on one or more of the other splitting categories in the composite split" << endl ;
	  delete obc ;
	  delete iter ;
	  return 0 ;
	}
      }
      delete iter ;

      // If a constrained split is specified, check that split parameter is a real-valued type
      if (pariter->second.second.size()>0) {
	if (!dynamic_cast<RooAbsReal*>(farg)) {
	  coutE(InputArguments) << "RooSimWSTool::build(" << GetName() << ") ERROR: Constrained split specified in non real-valued parameter " << farg->GetName() << endl ;
	  delete obc ;
	  return 0 ;
	}
      }

      // Fill object build config with object split rule
      osr._paramSplitMap[farg].first.add(splitCatSet) ;
      osr._paramSplitMap[farg].second = pariter->second.second ;

      // For multi-pdf configurations, check that the master index state name associated with this p.d.f exists as a state in the master category
      if (obc->_masterCat) {
	list<string>::iterator misi ;
	for (misi=sr._miStateNameList.begin() ; misi!=sr._miStateNameList.end() ; ++misi) {
	  const RooCatType* ctype = obc->_masterCat->lookupType(misi->c_str(),kFALSE) ;
	  if (ctype==0) {	  
	    coutE(ObjectHandling) << "RooSimWSTool::build(" << GetName() << ") ERROR: master index category " << obc->_masterCat->GetName() 
				  << " does not have a state named " << *misi << " which was specified as state associated with p.d.f " 
				  << sr.GetName() << endl ;
	    delete obc ;
	    return 0 ;	    
	  }
	  osr._miStateList.push_back(ctype) ;
	}	
      }

      // Add specified split cats to global list of all splitting categories
      obc->_usedSplitCats.add(splitCatSet,kTRUE) ;
      
    }
    // Need to add clause here for SplitRules without any split (which can happen in MultiBuildConfigs)
    if (sr._paramSplitMap.size()==0) {

      if (obc->_masterCat) {
	list<string>::iterator misi ;
	for (misi=sr._miStateNameList.begin() ; misi!=sr._miStateNameList.end() ; ++misi) {
	  const RooCatType* ctype = obc->_masterCat->lookupType(misi->c_str(),kFALSE) ;
	  if (ctype==0) {	  
	    coutE(ObjectHandling) << "RooSimWSTool::build(" << GetName() << ") ERROR: master index category " << obc->_masterCat->GetName() 
				  << " does not have a state named " << *misi << " which was specified as state associated with p.d.f " 
				  << sr.GetName() << endl ;
	    delete obc ;
	    return 0 ;	    
	  }
	  osr._miStateList.push_back(ctype) ;
	}	
      }
    }

    obc->_pdfmap[pdf] = osr ;

  }

  // Check validity of build restriction specifications, if any
  map<string,string>::iterator riter ;
  for (riter=bc._restr.begin() ; riter!=bc._restr.end() ; ++riter) {
    RooCategory* cat = _ws->cat(riter->first.c_str()) ;
    if (!cat) {
      coutE(ObjectHandling) << "RooSimWSTool::build(" << GetName() << ") ERROR: associated workspace " << _ws->GetName() 
			    << " does not contain a category named " << riter->first
			    << " for which build was requested to be restricted to states " << riter->second << endl ;
      delete obc ;
      return 0 ;
    }	

    char buf[4096] ;
    list<const RooCatType*> rlist ;
    strlcpy(buf,riter->second.c_str(),4096) ;
    
    char* tok = strtok(buf,"{,}") ;
    while(tok) {
      const RooCatType* ctype = cat->lookupType(tok,kFALSE) ;
      if (!ctype) {
	coutE(ObjectHandling) << "RooSimWSTool::build(" << GetName() << ") ERROR: restricted build category " << cat->GetName() 
			      << " does not have state " << tok << " as specified in restriction list" << endl ;	
	delete obc ;
	return 0 ;
      }
      rlist.push_back(ctype) ;
      tok = strtok(0,"{,}") ;
    }
    
    obc->_restr[cat] = rlist ;
  }
   
  return obc ;
}




////////////////////////////////////////////////////////////////////////////////
/// Internal build driver from validation ObjBuildConfig.

RooSimultaneous* RooSimWSTool::executeBuild(const char* simPdfName, ObjBuildConfig& obc, Bool_t verbose)
{
  RooArgSet cleanupList ;

  RooAbsCategoryLValue* physCat = obc._masterCat ;

  RooArgSet physModelSet ;
  map<string,RooAbsPdf*> stateMap ;

  map<RooAbsPdf*,ObjSplitRule>::iterator physIter = obc._pdfmap.begin() ;
  while(physIter!=obc._pdfmap.end()) {

    
    RooAbsPdf* physModel = physIter->first ;
    physModelSet.add(*physModel,kTRUE) ; // silence duplicate insertion warnings

    list<const RooCatType*>::iterator stiter ;
    for (stiter=physIter->second._miStateList.begin() ; stiter!=physIter->second._miStateList.end() ; ++stiter) {
      stateMap[(*stiter)->GetName()] = physModel ;
    }

    // Continue with next mapping
    ++physIter ;
  }
  if (verbose) coutI(ObjectHandling) << "RooSimWSTool::executeBuild: list of prototype pdfs " << physModelSet << endl ;

  RooArgSet splitCatSet(obc._usedSplitCats) ;
  if (physCat) splitCatSet.add(*physCat) ;

  RooArgSet splitCatSetFund ;
  TIterator* scsiter = splitCatSet.createIterator() ;
  RooAbsCategory* scat ;
  while((scat=(RooAbsCategory*)scsiter->Next())) {
    if (scat->isFundamental()) {
      splitCatSetFund.add(*scat) ;
    } else {
      RooArgSet* scatvars = scat->getVariables() ;
      splitCatSetFund.add(*scatvars) ;
      delete scatvars ;
    }
  }
  delete scsiter ;


  RooAbsCategoryLValue* masterSplitCat ;
  if (splitCatSetFund.getSize()>1) {
    masterSplitCat = new RooSuperCategory("masterSplitCat","Master splitting category",splitCatSetFund) ;
  } else {
    masterSplitCat = (RooAbsCategoryLValue*) splitCatSetFund.first() ;
  }
  if (verbose) coutI(ObjectHandling) << "RooSimWSTool::executeBuild: list of splitting categories " << splitCatSet << endl ;

  RooArgSet splitNodeListOwned ; // owns all newly created components
  RooArgSet splitNodeListAll ; // all leaf nodes, preload with ws contents to auto-connect existing specializations
  TList* customizerList = new TList ;

  // Loop over requested physics models and build components
  TIterator* physMIter = physModelSet.createIterator() ;
  RooAbsPdf* physModel ;
  while((physModel=(RooAbsPdf*)physMIter->Next())) {
    if (verbose) coutI(ObjectHandling) << "RooSimPdfBuilder::executeBuild: processing prototype pdf " << physModel->GetName() << endl ;

    RooCustomizer* physCustomizer = new RooCustomizer(*physModel,*masterSplitCat,splitNodeListOwned,&splitNodeListAll) ;
    customizerList->Add(physCustomizer) ;

    map<RooAbsArg*, pair<RooArgSet,string> >::iterator splitIter ;
    for (splitIter = obc._pdfmap[physModel]._paramSplitMap.begin() ; splitIter != obc._pdfmap[physModel]._paramSplitMap.end() ; ++splitIter) {

      // If split is composite, first make multicategory with name 'A,B,C' and insert in WS
      
      // Construct name of (composite) split category (function)
      RooArgSet& splitCatSetTmp = splitIter->second.first ;
      string splitName = makeSplitName(splitCatSetTmp) ;

      // If composite split object does not exist yet, create it now
      RooAbsCategory* splitCat = _ws->catfunc(splitName.c_str()) ;
      if (!splitCat) {
	splitCat = new RooMultiCategory(splitName.c_str(),splitName.c_str(),splitCatSetTmp) ;
	cleanupList.addOwned(*splitCat) ;
	_ws->import(*splitCat,RooFit::Silence(!verbose)) ;
      }
            
      // If remainder category needs to be made, create RFV of appropriate for that and insert in WS
      if(splitIter->second.second.size()>0) {
	
	// Check that specified split name is in fact valid
	if (!splitCat->hasLabel(splitIter->second.second)) {
	  coutE(InputArguments) << "RooSimWSTool::executeBuild(" << GetName() << ") ERROR: name of remainder state for constrained split, '" 
				<< splitIter->second.second << "' , does not match any state name of (composite) split category " << splitCat->GetName() << endl ;
	  return 0 ;
	}

	// First build manually the specializations of all non-remainder states, as the remainder state depends on these
	RooArgSet fracLeafList ;
	TIterator* sctiter = splitCat->typeIterator() ;
	RooCatType* type ;
	while((type=(RooCatType*)sctiter->Next())) {
	  
	  // Skip remainder state
	  if (splitIter->second.second == type->GetName()) continue ;
	  	  
	  // Construct name of split leaf
	  TString splitLeafName(splitIter->first->GetName()) ;
	  splitLeafName.Append("_") ;
	  splitLeafName.Append(type->GetName()) ;
	  
	  // Check if split leaf already exists	  
	  RooAbsArg* splitLeaf = _ws->fundArg(splitLeafName) ;
	  if (!splitLeaf) {
	    // If not create it now
	    splitLeaf = (RooAbsArg*) splitIter->first->clone(splitLeafName) ;
	    _ws->import(*splitLeaf,RooFit::Silence(!verbose)) ;
	  }
	  fracLeafList.add(*splitLeaf) ;
	}
	delete sctiter ;		
	

	// Build specialization for remainder state and insert in workspace 
	RooFracRemainder* fracRem = new RooFracRemainder(Form("%s_%s",splitIter->first->GetName(),splitIter->second.second.c_str()),"Remainder fraction",fracLeafList) ;
	cleanupList.addOwned(*fracRem) ;
	_ws->import(*fracRem) ;

      }
      

      // Add split definition to customizer
      physCustomizer->splitArgs(*splitIter->first,*splitCat) ;
    }
  }
  delete physMIter ;

  // List all existing workspace components as prebuilt items for the customizers at this point
  splitNodeListAll.add(_ws->components()) ;

  if (verbose) coutI(ObjectHandling)  << "RooSimWSTool::executeBuild: configured customizers for all prototype pdfs" << endl ;

  // Create fit category from physCat and splitCatList ;
  RooArgSet fitCatList ;
  if (physCat) fitCatList.add(*physCat) ;

  // Add observables of splitCatSet members, rather than splitCatSet members directly
  // as there may be cat->cat functions in here
  scsiter = splitCatSet.createIterator() ;
  while((scat=(RooAbsCategory*)scsiter->Next())) {
    if (scat->isFundamental()) {
      fitCatList.add(*scat) ;
    } else {
      RooArgSet* scatvars = scat->getVariables() ;
      fitCatList.add(*scatvars) ;
      delete scatvars ;
    }
  }
  delete scsiter ;


  TIterator* fclIter = fitCatList.createIterator() ;
  string mcatname = string(simPdfName) + "_index" ;
  RooAbsCategoryLValue* fitCat = 0 ;
  if (fitCatList.getSize()>1) {
    fitCat = new RooSuperCategory(mcatname.c_str(),mcatname.c_str(),fitCatList) ;
    cleanupList.addOwned(*fitCat) ;
  } else {
    fitCat = (RooAbsCategoryLValue*) fitCatList.first() ;
  }

  // Create master PDF 
  RooSimultaneous* simPdf = new RooSimultaneous(simPdfName,simPdfName,*fitCat) ;
  cleanupList.addOwned(*simPdf) ;
  
  // Add component PDFs to master PDF
  TIterator* fcIter = fitCat->typeIterator() ;

  RooCatType* fcState ;  
  while((fcState=(RooCatType*)fcIter->Next())) {
    // Select fitCat state
    fitCat->setLabel(fcState->GetName()) ;

    // Check if this fitCat state is selected
    fclIter->Reset() ;
    RooAbsCategory* splitCat ;
    Bool_t select(kFALSE) ;
    if (obc._restr.size()>0) {
      while((splitCat=(RooAbsCategory*)fclIter->Next())) {
	// Find selected state list 
	
	list<const RooCatType*> slist = obc._restr[splitCat] ;    
	if (slist.empty()) {
	  continue ;
	}
	
	list<const RooCatType*>::iterator sli ;
	for (sli=slist.begin() ; sli!=slist.end() ; ++sli) {
	  if (string(splitCat->getCurrentLabel())==(*sli)->GetName()) {
	    select=kTRUE ;
	  }
	}      
      }
      if (!select) continue ;
    } else {
      select = kTRUE ;
    }
    
    // Select appropriate PDF for this physCat state
    RooCustomizer* physCustomizer ;
    if (physCat) {      
      RooAbsPdf* pdf = stateMap[physCat->getCurrentLabel()] ;
      if (pdf==0) {
	continue ;
      }
      physCustomizer = (RooCustomizer*) customizerList->FindObject(pdf->GetName());  
    } else {
      physCustomizer = (RooCustomizer*) customizerList->First() ;
    }

    if (verbose) coutI(ObjectHandling) << "RooSimWSTool::executeBuild: Customizing prototype pdf " << physCustomizer->GetName() 
				       << " for mode " << fcState->GetName() << endl ;    

    // Customizer PDF for current state and add to master simPdf
    RooAbsPdf* fcPdf = (RooAbsPdf*) physCustomizer->build(masterSplitCat->getCurrentLabel(),kFALSE) ;
    simPdf->addPdf(*fcPdf,fcState->GetName()) ;
  }
  delete fcIter ;
  
  _ws->import(*simPdf,obc._conflProtocol,RooFit::Silence(!verbose)) ;

  // Delete customizers
  customizerList->Delete() ;
  delete customizerList ;
  delete fclIter ;
  return (RooSimultaneous*) _ws->pdf(simPdf->GetName()) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Construct name of composite split

std::string RooSimWSTool::makeSplitName(const RooArgSet& splitCatSet) 
{
  string name ;

  TIterator* iter = splitCatSet.createIterator() ;
  RooAbsArg* arg ;
  Bool_t first=kTRUE ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (first) {
      first=kFALSE;
    } else {
      name += "," ;
    }
    name += arg->GetName() ;
  }
  delete iter ;

  return name ;
}




////////////////////////////////////////////////////////////////////////////////
/// Specify that parameters names listed in paramNameList be split in (product of) category(s)
/// listed in categoryNameList

void RooSimWSTool::SplitRule::splitParameter(const char* paramNameList, const char* categoryNameList) 
{
  char paramBuf[4096] ;
  char catBuf[4096] ;
  strlcpy(paramBuf,paramNameList,4096) ;
  strlcpy(catBuf,categoryNameList,4096) ;

  // First parse category list
  list<string> catList ;
  char* cat = strtok(catBuf,"{,}") ;
  while(cat) {
    catList.push_back(cat) ;
    cat = strtok(0,"{,}") ;
  }

  // Now parse parameter list
  char* param = strtok(paramBuf,"{,}") ;
  while(param) {
    _paramSplitMap[param] = pair<list<string>,string>(catList,"") ;
    param = strtok(0,"{,}") ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Specify that parameters names listed in paramNameList be split in constrained way in (product of) category(s)
/// listed in categoryNameList and that remainder fraction formula be put in state with name remainderStateName

void RooSimWSTool::SplitRule::splitParameterConstrained(const char* paramNameList, const char* categoryNameList, const char* remainderStateName) 
{
  char paramBuf[4096] ;
  char catBuf[4096] ;
  strlcpy(paramBuf,paramNameList,4096) ;
  strlcpy(catBuf,categoryNameList,4096) ;

  // First parse category list
  list<string> catList ;
  char* cat = strtok(catBuf,"{,}") ;
  while(cat) {
    catList.push_back(cat) ;
    cat = strtok(0,"{,}") ;
  }

  // Now parse parameter list
  char* param = strtok(paramBuf,"{,}") ;
  while(param) {
    _paramSplitMap[param] = pair<list<string>,string>(catList,remainderStateName) ;
    param = strtok(0,"{,}") ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Construct the SplitRule object from a list of named arguments past to RooSimWSTool::build
/// This method parses any SplitParam and SplitParamComstrained argument in the list

void RooSimWSTool::SplitRule::configure(const RooCmdArg& arg1,const RooCmdArg& arg2,const RooCmdArg& arg3,
					   const RooCmdArg& arg4, const RooCmdArg& arg5,const RooCmdArg& arg6)
{
  list<const RooCmdArg*> cmdList ;  
  cmdList.push_back(&arg1) ;  cmdList.push_back(&arg2) ;
  cmdList.push_back(&arg3) ;  cmdList.push_back(&arg4) ;
  cmdList.push_back(&arg5) ;  cmdList.push_back(&arg6) ;

  list<const RooCmdArg*>::iterator iter ;
  for (iter=cmdList.begin() ; iter!=cmdList.end() ; ++iter) {

    if ((*iter)->opcode()==0) continue ;

    string name = (*iter)->opcode() ;

    if (name=="SplitParam") {
      splitParameter((*iter)->getString(0),(*iter)->getString(1)) ;
    } else if (name=="SplitParamConstrained") {
      splitParameterConstrained((*iter)->getString(0),(*iter)->getString(1),(*iter)->getString(2)) ;
    }
  }
}




////////////////////////////////////////////////////////////////////////////////
/// Add prototype p.d.f pdfName to build configuration with associated split rules 'sr'

RooSimWSTool::BuildConfig::BuildConfig(const char* pdfName, SplitRule& sr)
{
  internalAddPdf(pdfName,"",sr) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Construct build configuration from single prototype 'pdfName' and list of arguments
/// that can be passed to RooSimWSTool::build() method. This routine parses SplitParam()
/// SplitParamConstrained() and Restrict() arguments.

RooSimWSTool::BuildConfig::BuildConfig(const char* pdfName, const RooCmdArg& arg1,const RooCmdArg& arg2,
					  const RooCmdArg& arg3,const RooCmdArg& arg4, const RooCmdArg& arg5,const RooCmdArg& arg6)
{
  SplitRule sr(pdfName) ;
  sr.configure(arg1,arg2,arg3,arg4,arg5,arg6) ;
  internalAddPdf(pdfName,"",sr) ;
  _conflProtocol = RooFit::RenameConflictNodes(pdfName) ;

  list<const RooCmdArg*> cmdList ;  
  cmdList.push_back(&arg1) ;  cmdList.push_back(&arg2) ;
  cmdList.push_back(&arg3) ;  cmdList.push_back(&arg4) ;
  cmdList.push_back(&arg5) ;  cmdList.push_back(&arg6) ;

  list<const RooCmdArg*>::iterator iter ;
  for (iter=cmdList.begin() ; iter!=cmdList.end() ; ++iter) {
    if ((*iter)->opcode()==0) continue ;
    string name = (*iter)->opcode() ;
    if (name=="Restrict") {
      restrictBuild((*iter)->getString(0),(*iter)->getString(1)) ;
    }    
    if (name=="RenameConflictNodes") {
      _conflProtocol = *(*iter) ;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor to make BuildConfig from legacy RooSimPdfBuilder configuration
/// Empty for now

RooSimWSTool::BuildConfig::BuildConfig(const RooArgSet& /*legacyBuildConfig*/) 
{
}


////////////////////////////////////////////////////////////////////////////////
/// Internal routine to add prototype pdf 'pdfName' with list of associated master states 'miStateNameList
/// and split rules 'sr' to configuration

void RooSimWSTool::BuildConfig::internalAddPdf(const char* pdfName, const char* miStateNameList,SplitRule& sr) 
{
  char buf[4096] ;
  strlcpy(buf,miStateNameList,4096) ;
  
  char* tok = strtok(buf,",") ;
  while(tok) {
    sr._miStateNameList.push_back(tok) ;
    tok = strtok(0,",") ;
  }
    
  _pdfmap[pdfName] = sr ;  
}


////////////////////////////////////////////////////////////////////////////////
/// Restrict build by only considering state names in stateList for split in category catName

void RooSimWSTool::BuildConfig::restrictBuild(const char* catName, const char* stateList) 
{
  _restr[catName] = stateList ;
}




////////////////////////////////////////////////////////////////////////////////
/// Construct MultiBuildConfig for build configuration with multiple prototype p.d.f.s
/// masterIndexCat is the name of the master index category that decides which
/// prototype is used.

RooSimWSTool::MultiBuildConfig::MultiBuildConfig(const char* masterIndexCat)  
{
  _masterCatName = masterIndexCat ;
}



////////////////////////////////////////////////////////////////////////////////
/// Add protytpe p.d.f 'pdfName' to MultiBuildConfig associated with master indes states 'miStateList'. This
/// method parses the SplitParam() and SplitParamConstrained() arguments

void RooSimWSTool::MultiBuildConfig::addPdf(const char* miStateList, const char* pdfName, const RooCmdArg& arg1,const RooCmdArg& arg2,
					       const RooCmdArg& arg3,const RooCmdArg& arg4, const RooCmdArg& arg5,const RooCmdArg& arg6)
{
  SplitRule sr(pdfName) ;
  sr.configure(arg1,arg2,arg3,arg4,arg5,arg6) ;
  internalAddPdf(pdfName,miStateList,sr) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Add protytpe p.d.f 'pdfName' to MultiBuildConfig associated with master indes states 'miStateList'. 

void RooSimWSTool::MultiBuildConfig::addPdf(const char* miStateList, const char* pdfName, SplitRule& sr) 
{
  internalAddPdf(pdfName,miStateList,sr) ;
}




////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooSimWSTool::ObjSplitRule::~ObjSplitRule()
{
}




////////////////////////////////////////////////////////////////////////////////
/// Print details of a validated build configuration

void RooSimWSTool::ObjBuildConfig::print()
{
  // --- Dump contents of object build config ---
  map<RooAbsPdf*,ObjSplitRule>::iterator ri ;
  for (ri = _pdfmap.begin() ; ri != _pdfmap.end() ; ++ri ) {    
    cout << "Splitrule for p.d.f " << ri->first->GetName() << " with state list " ;
    for (std::list<const RooCatType*>::iterator misi= ri->second._miStateList.begin() ; misi!=ri->second._miStateList.end() ; ++misi) {
      cout << (*misi)->GetName() << " " ;
    }
    cout << endl ;

    map<RooAbsArg*,pair<RooArgSet,string> >::iterator csi ;
    for (csi = ri->second._paramSplitMap.begin() ; csi != ri->second._paramSplitMap.end() ; ++csi ) {    
      if (csi->second.second.length()>0) {
	cout << " parameter " << csi->first->GetName() << " is split with constraint in categories " << csi->second.first 
	     << " with remainder in state " << csi->second.second << endl ;      
      } else {
	cout << " parameter " << csi->first->GetName() << " is split with constraint in categories " << csi->second.first << endl ;      
      }
    }        
  }

  map<RooAbsCategory*,list<const RooCatType*> >::iterator riter ;
  for (riter=_restr.begin() ; riter!=_restr.end() ; ++riter) {
    cout << "Restricting build in category " << riter->first->GetName() << " to states " ;
    list<const RooCatType*>::iterator i ;
    for (i=riter->second.begin() ; i!=riter->second.end() ; ++i) {
      if (i!=riter->second.begin()) cout << "," ;
      cout << (*i)->GetName() ;
    }
    cout << endl ;
  }

}




////////////////////////////////////////////////////////////////////////////////

std::string RooSimWSTool::SimWSIFace::create(RooFactoryWSTool& ft, const char* typeName, const char* instanceName, std::vector<std::string> args) 
{
  string tn(typeName) ;
  if (tn=="SIMCLONE") {

    // Perform syntax check. Warn about any meta parameters other than $SplitParam, $SplitParamConstrained, $Restrict and $Verbose
    for (unsigned int i=1 ; i<args.size() ; i++) {
      if (args[i].find("$SplitParam(")!=0 &&
	  args[i].find("$SplitParamConstrained(")!=0 &&
	  args[i].find("$SplitRestrict(")!=0 &&
	  args[i].find("$Verbose(")!=0) {
	throw string(Form("RooSimWSTool::SimWSIFace::create() ERROR: unknown token %s encountered",args[i].c_str())) ;
      }
    }

    // Make SplitRule object from $SplitParam and $SplitParamConstrained arguments
    RooSimWSTool::SplitRule sr(args[0].c_str()) ;
    for (unsigned int i=1 ; i<args.size() ; i++) {
      if (args[i].find("$SplitParam(")==0) {
	vector<string> subargs = ft.splitFunctionArgs(args[i].c_str()) ;
	if (subargs.size()!=2) {
	  throw string(Form("Incorrect number of arguments in $SplitParam, have %d, expect 2",(Int_t)subargs.size())) ;
	}
	sr.splitParameter(subargs[0].c_str(),subargs[1].c_str()) ;
      } else if (args[i].find("$SplitParamConstrained(")==0) {
	vector<string> subargs = ft.splitFunctionArgs(args[i].c_str()) ;
	if (subargs.size()!=3) {
	  throw string(Form("Incorrect number of arguments in $SplitParamConstrained, have %d, expect 3",(Int_t)subargs.size())) ;
	}
	sr.splitParameterConstrained(subargs[0].c_str(), subargs[1].c_str(), subargs[2].c_str()) ;	
      } 
    }

    // Make BuildConfig object
    RooSimWSTool::BuildConfig bc(args[0].c_str(),sr) ;
    for (unsigned int i=1 ; i<args.size() ; i++) {
      if (args[i].find("$Restrict(")==0) {
	vector<string> subargs = ft.splitFunctionArgs(args[i].c_str()) ;
	if (subargs.size()!=2) {
	  throw string(Form("Incorrect number of arguments in $Restrict, have %d, expect 2",(Int_t)subargs.size())) ;
	}
	bc.restrictBuild(subargs[0].c_str(),subargs[1].c_str()) ;
      }
    }

    // Look for verbose flag
    Bool_t verbose(kFALSE) ;
    for (unsigned int i=1 ; i<args.size() ; i++) {
      if (args[i].find("$Verbose(")==0) {
	vector<string> subargs = ft.splitFunctionArgs(args[i].c_str()) ;
	if (subargs.size()>0) {
	  verbose = atoi(subargs[0].c_str()) ;
	}
      }
    }

    // Build pdf clone
    RooSimWSTool sct(ft.ws()) ;    
    RooAbsPdf* pdf = sct.build(instanceName,bc,verbose) ;
    if (!pdf) {
      throw string(Form("RooSimWSTool::SimWSIFace::create() error in RooSimWSTool::build() for %s",instanceName)) ;
    }
    
    // Import into workspace
    ft.ws().import(*pdf,RooFit::Silence()) ;

  } else if (tn=="MSIMCLONE") {

    // First make a multibuild config from the master index cat
    RooSimWSTool::MultiBuildConfig mbc(args[0].c_str()) ;

    for (unsigned int i=1 ; i<args.size() ; i++) {
      if (args[i].find("$AddPdf(")==0) {
	// Process an add-pdf operation
	vector<string> subargs = ft.splitFunctionArgs(args[i].c_str()) ;

	// Make SplitRule object from $SplitParam and $SplitParamConstrained arguments
	RooSimWSTool::SplitRule sr(subargs[1].c_str()) ;
	for (unsigned int j=2 ; j<subargs.size() ; j++) {
	  if (subargs[j].find("$SplitParam(")==0) {
	    vector<string> subsubargs = ft.splitFunctionArgs(subargs[j].c_str()) ;
	    if (subsubargs.size()!=2) {
	      throw string(Form("Incorrect number of arguments in $SplitParam, have %d, expect 2",(Int_t)subsubargs.size())) ;
	    }
	    sr.splitParameter(subsubargs[0].c_str(),subsubargs[1].c_str()) ;
	  } else if (subargs[j].find("$SplitParamConstrained(")==0) {
	    vector<string> subsubargs = ft.splitFunctionArgs(subargs[j].c_str()) ;
	    if (subsubargs.size()!=3) {
	      throw string(Form("Incorrect number of arguments in $SplitParamConstrained, have %d, expect 3",(Int_t)subsubargs.size())) ;
	    }
	    sr.splitParameterConstrained(subsubargs[0].c_str(), subsubargs[1].c_str(), subsubargs[2].c_str()) ;	
	  } 
	}     
	mbc.addPdf(subargs[0].c_str(),subargs[1].c_str(),sr) ;

      } else if (args[i].find("$Restrict(")==0) {
	
	// Process a restrict operation 	
	vector<string> subargs = ft.splitFunctionArgs(args[i].c_str()) ;
	if (subargs.size()!=2) {
	  throw string(Form("Incorrect number of arguments in $Restrict, have %d, expect 2",(Int_t)subargs.size())) ;
	}
	mbc.restrictBuild(subargs[0].c_str(),subargs[1].c_str()) ;
	
      } else {
	throw string(Form("RooSimWSTool::SimWSIFace::create() ERROR: unknown token in MSIMCLONE: %s",args[i].c_str())) ;
      }
    }    

    // Build pdf clone
    RooSimWSTool sct(ft.ws()) ;    
    RooAbsPdf* pdf = sct.build(instanceName,mbc,kFALSE) ;
    if (!pdf) {
      throw string(Form("RooSimWSTool::SimWSIFace::create() error in RooSimWSTool::build() for %s",instanceName)) ;
    }
    
    // Import into workspace
    ft.ws().import(*pdf,RooFit::Silence()) ;

    
  } else {
    throw string(Form("RooSimWSTool::SimWSIFace::create() ERROR: Unknown meta-type %s requested",typeName)) ;
  }
  
  return string(instanceName) ;
}
