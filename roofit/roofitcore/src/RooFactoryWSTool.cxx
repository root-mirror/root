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
\file RooFactoryWSTool.cxx
\class RooFactoryWSTool
\ingroup Roofitcore

RooFactoryWSTool is a class similar to TTree::MakeClass() that generates
skeleton code for RooAbsPdf and RooAbsReal functions given
a list of input parameter names. The factory can also compile
the generated code on the fly, and on request also
instantiate the objects.

It interprets all expressions for RooWorkspace::factory(const char*).
**/

#include "RooFit.h"

#include "RooFactoryWSTool.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooMsgService.h"
#include "RooWorkspace.h"
#include "TInterpreter.h"
#include "TEnum.h"
#include "RooAbsPdf.h"
#include <fstream>
#include "strtok.h"
#include "strlcpy.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooSimultaneous.h"
#include "RooFFTConvPdf.h"
#include "RooNumConvPdf.h"
#include "RooResolutionModel.h"
#include "RooProduct.h"
#include "RooAddition.h"
#include "RooChi2Var.h"
#include "RooNLLVar.h"
#include "RooRealSumPdf.h"
#include "RooConstVar.h"
#include "RooDerivative.h"
#include "RooStringVar.h"
#include "TROOT.h"

using namespace RooFit ;
using namespace std ;

#define BUFFER_SIZE 64000 

ClassImp(RooFactoryWSTool); 
;

RooFactoryWSTool* RooFactoryWSTool::_of = 0 ;
map<string,RooFactoryWSTool::IFace*>* RooFactoryWSTool::_hooks=0 ;

namespace {

static Int_t init();

Int_t dummy = init() ;

static Int_t init()
{
  RooFactoryWSTool::IFace* iface = new RooFactoryWSTool::SpecialsIFace ;

  // Operator p.d.f.s
  RooFactoryWSTool::registerSpecial("SUM",iface) ;
  RooFactoryWSTool::registerSpecial("RSUM",iface) ;
  RooFactoryWSTool::registerSpecial("ASUM",iface) ;
  RooFactoryWSTool::registerSpecial("PROD",iface) ;
  RooFactoryWSTool::registerSpecial("SIMUL",iface) ;
  RooFactoryWSTool::registerSpecial("EXPR",iface) ;
  RooFactoryWSTool::registerSpecial("FCONV",iface) ;
  RooFactoryWSTool::registerSpecial("NCONV",iface) ;

  // Operator functions
  RooFactoryWSTool::registerSpecial("sum",iface) ;
  RooFactoryWSTool::registerSpecial("prod",iface) ;
  RooFactoryWSTool::registerSpecial("expr",iface) ;
  RooFactoryWSTool::registerSpecial("nconv",iface) ;

  // Test statistics
  RooFactoryWSTool::registerSpecial("nll",iface) ;
  RooFactoryWSTool::registerSpecial("chi2",iface) ;
  RooFactoryWSTool::registerSpecial("profile",iface) ;

  // Integration and derivation
  RooFactoryWSTool::registerSpecial("int",iface) ;
  RooFactoryWSTool::registerSpecial("deriv",iface) ;
  RooFactoryWSTool::registerSpecial("cdf",iface) ;
  RooFactoryWSTool::registerSpecial("PROJ",iface) ;

  // Miscellaneous
  RooFactoryWSTool::registerSpecial("dataobs",iface) ;
  RooFactoryWSTool::registerSpecial("set",iface) ;

  (void) dummy;
  return 0 ;
}

}

#ifndef _WIN32
#include <strings.h>
#endif



////////////////////////////////////////////////////////////////////////////////

RooFactoryWSTool::RooFactoryWSTool(RooWorkspace& inws) : _ws(&inws), _errorCount(0), _autoClassPostFix("")

{
  // Default constructor  
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooFactoryWSTool::~RooFactoryWSTool() 
{
}




////////////////////////////////////////////////////////////////////////////////
/// Low-level factory interface for creating a RooRealVar with a given range and initial value

RooRealVar* RooFactoryWSTool::createVariable(const char* name, Double_t xmin, Double_t xmax)
{
  // First check if variable already exists
  if (_ws->var(name)) {
    coutE(ObjectHandling) << "RooFactoryWSTool::createFactory() ERROR: variable with name '" << name << "' already exists" << endl ;
    logError() ;
    return 0 ;
  }

  // Create variable
  RooRealVar var(name,name,xmin,xmax) ;

  // Put in workspace
  if (_ws->import(var,Silence())) logError() ;

  return _ws->var(name) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Low-level factory interface for creating a RooCategory with a given list of state names. The State name list
/// can be of the form 'name1,name2,name3' or of the form 'name1=id1,name2=id2,name3=id3'

RooCategory* RooFactoryWSTool::createCategory(const char* name, const char* stateNameList) 
{
  // Create variable
  RooCategory cat(name,name) ;

  // Add listed state names
  if (stateNameList) {
     const size_t tmpSize = strlen(stateNameList)+1;
    char *tmp = new char[tmpSize] ;
    strlcpy(tmp,stateNameList,tmpSize) ;
    char* save ;
    char* tok = R__STRTOK_R(tmp,",",&save) ;
    while(tok) {
      char* sep = strchr(tok,'=') ;
      if (sep) {
	*sep = 0 ;
	Int_t id = atoi(sep+1) ;
	cat.defineType(tok,id) ;
	*sep = '=' ;
      } else {
	cat.defineType(tok) ;
      }
      tok = R__STRTOK_R(0,",",&save) ;
    }
    delete[] tmp ;
  }

  cat.setStringAttribute("factory_tag",Form("%s[%s]",name,stateNameList)) ;

  // Put in workspace
  if (_ws->import(cat,Silence())) logError() ;  

  return _ws->cat(name) ;
}

namespace {
  static bool isEnum(const char* classname) {
    // Returns true if given type is an enum
    ClassInfo_t* cls = gInterpreter->ClassInfo_Factory(classname);
    long property = gInterpreter->ClassInfo_Property(cls);
    gInterpreter->ClassInfo_Delete(cls);
    return (property&kIsEnum);
  }


  static bool isValidEnumValue(const char* enumName, const char* enumConstantName) {
    // Returns true if given type is an enum

    if (!enumName) return false;

    auto theEnum = TEnum::GetEnum(enumName);
    if (!enumName) return false;

    // Attempt 1: Enum constant name as is
    if (theEnum->GetConstant(enumConstantName)) return true;
    // Attempt 2: Remove the scope preceding the enum constant name
    auto tmp = strstr(enumConstantName, "::");
    if (tmp) {
      auto enumConstantNameNoScope = tmp+2;
      if (theEnum->GetConstant(enumConstantNameNoScope)) return true;
    }

    return false;
  }

  static pair<list<string>,unsigned int> ctorArgs(const char* classname, UInt_t nMinArg) {
    // Utility function for RooFactoryWSTool. Return arguments of 'first' non-default, non-copy constructor of any RooAbsArg
    // derived class. Only constructors that start with two 'const char*' arguments (for name and title) are considered
    // The returned object contains

    Int_t nreq(0);
    list<string> ret;

    ClassInfo_t* cls = gInterpreter->ClassInfo_Factory(classname);
    MethodInfo_t* func = gInterpreter->MethodInfo_Factory(cls);
    while(gInterpreter->MethodInfo_Next(func)) {
      ret.clear();
      nreq=0;

      // Find 'the' constructor

      // Skip non-public methods
      if (!(gInterpreter->MethodInfo_Property(func) & kIsPublic)) {
        continue;
      }

      // Return type must be class name
      if (string(classname) != gInterpreter->MethodInfo_TypeName(func)) {
        continue;
      }

      // Skip default constructor
      int nargs = gInterpreter->MethodInfo_NArg(func);
      if (nargs==0 || nargs==gInterpreter->MethodInfo_NDefaultArg(func)) {
        continue;
      }

      MethodArgInfo_t* arg = gInterpreter->MethodArgInfo_Factory(func);
      while (gInterpreter->MethodArgInfo_Next(arg)) {
        // Require that first two arguments are of type const char*
        const char* argTypeName = gInterpreter->MethodArgInfo_TypeName(arg);
        if (nreq<2 && ((string("char*") != argTypeName
                && !(gInterpreter->MethodArgInfo_Property(arg) & kIsConstPointer))
              && string("const char*") != argTypeName)) {
          continue ;
        }
        ret.push_back(argTypeName) ;
        if(!gInterpreter->MethodArgInfo_DefaultValue(arg)) nreq++;
      }
      gInterpreter->MethodArgInfo_Delete(arg);

      // Check that the number of required arguments is at least nMinArg
      if (ret.size()<nMinArg) {
        continue;
      }

      break;
    }
    gInterpreter->MethodInfo_Delete(func);
    gInterpreter->ClassInfo_Delete(cls);
    return pair<list<string>,unsigned int>(ret,nreq);
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Low-level factory interface for creating a RooAbsPdf of a given class with a given list of input variables
/// The variable list varList should be of the form "a,b,c" where the interpretation of the argument is
/// dependent on the p.d.f. Set and List arguments can be passed by substituting a single argument with
/// the form (a,b,c), i.e. one can set varList to "x,(a0,a1,a2)" to pass a RooAbsReal and a RooArgSet as arguments.

RooAbsArg* RooFactoryWSTool::createArg(const char* className, const char* objName, const char* varList) 
{
  // Find class in ROOT class table
  TClass* tc = resolveClassName(className);
  if (!tc) {
    coutE(ObjectHandling) << "RooFactoryWSTool::createArg() ERROR class " << className << " not found in factory alias table, nor in ROOT class table" << endl;
    logError();
    return 0;
  }

  className = tc->GetName();

  // Check that class inherits from RooAbsPdf
  if (!tc->InheritsFrom(RooAbsArg::Class())) {
    coutE(ObjectHandling) << "RooFactoryWSTool::createArg() ERROR class " << className << " does not inherit from RooAbsArg" << endl;
    logError();
    return 0;
  }

  _args.clear();
  string tmp(varList);
  size_t blevel = 0, end_tok, start_tok = 0;
  bool litmode = false;
  for (end_tok = 0; end_tok < tmp.length(); end_tok++) {
    // Keep track of opening and closing brackets
    if (tmp[end_tok]=='{' || tmp[end_tok]=='(' || tmp[end_tok]=='[') blevel++;
    if (tmp[end_tok]=='}' || tmp[end_tok]==')' || tmp[end_tok]==']') blevel--;

    // Keep track of string literals
    if (tmp[end_tok]=='"' || tmp[end_tok]=='\'') litmode = !litmode;

    // If we encounter a comma at zero bracket level
    // push the current substring from start_tok to end_tok
    // and start the next token
    if (litmode == false && blevel == 0 && tmp[end_tok] == ',') {
      _args.push_back(tmp.substr(start_tok, end_tok - start_tok));
      start_tok = end_tok+1;
    }
  }
  _args.push_back(tmp.substr(start_tok, end_tok));

  // Try CINT interface
  pair<list<string>,unsigned int> ca = ctorArgs(className,_args.size()+2) ;
  if (ca.first.size()==0) {
    coutE(ObjectHandling) << "RooFactoryWSTool::createArg() ERROR no suitable constructor found for class " << className << endl ;
    logError() ;
    return 0 ;
  }  


  // Check if number of provided args is in valid range (add two to accomodate name and title strings)
  if (_args.size()+2<ca.second || _args.size()+2>ca.first.size()) {
    if (ca.second==ca.first.size()) {
      coutE(ObjectHandling) << "RooFactoryWSTool::createArg() ERROR number of arguments provided (" << _args.size() << ") for class is invalid, " << className 
			    << " expects " << ca.first.size()-2 << endl ;
      logError() ;
    } else {
      coutE(ObjectHandling) << "RooFactoryWSTool::createArg() ERROR number of arguments provided (" << _args.size() << ") for class is invalid " << className 
			    << " expect number between " << ca.second-2 << " and " << ca.first.size()-2 << endl ;
      logError() ;
    }
    return 0 ;
  }

  // Now construct CINT constructor spec, start with mandatory name and title args
  string cintExpr(Form("new %s(\"%s\",\"%s\"",className,objName,objName)) ;

  // Install argument in static data member to be accessed below through static CINT interface functions  
  _of = this ;

  
  try {
    Int_t i(0) ;
    list<string>::iterator ti = ca.first.begin() ; ++ti ; ++ti ;
    for (vector<string>::iterator ai = _args.begin() ; ai != _args.end() ; ++ai,++ti,++i) {
      if ((*ti)=="RooAbsReal&" || (*ti)=="const RooAbsReal&") {
	RooFactoryWSTool::as_FUNC(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_FUNC(%d)",i) ;
      } else if ((*ti)=="RooAbsArg&" || (*ti)=="const RooAbsArg&") {
	RooFactoryWSTool::as_ARG(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_ARG(%d)",i) ;
      } else if ((*ti)=="RooRealVar&" || (*ti)=="const RooRealVar&") {
	RooFactoryWSTool::as_VAR(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_VAR(%d)",i) ;
      } else if ((*ti)=="RooAbsRealLValue&" || (*ti)=="const RooAbsRealLValue&") {
	RooFactoryWSTool::as_VARLV(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_VARLV(%d)",i) ;
      } else if ((*ti)=="RooCategory&" || (*ti)=="const RooCategory&") {
	RooFactoryWSTool::as_CAT(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_CAT(%d)",i) ;
      } else if ((*ti)=="RooAbsCategory&" || (*ti)=="const RooAbsCategory&") {
	RooFactoryWSTool::as_CATFUNC(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_CATFUNC(%d)",i) ;
      } else if ((*ti)=="RooAbsCategoryLValue&" || (*ti)=="const RooAbsCategoryLValue&") {
	RooFactoryWSTool::as_CATLV(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_CATLV(%d)",i) ;
      } else if ((*ti)=="RooAbsPdf&" || (*ti)=="const RooAbsPdf&") {
	RooFactoryWSTool::as_PDF(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_PDF(%d)",i) ;
      } else if ((*ti)=="RooResolutionModel&" || (*ti)=="const RooResolutionModel&") {
	RooFactoryWSTool::as_RMODEL(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_RMODEL(%d)",i) ;
      } else if ((*ti)=="RooAbsData&" || (*ti)=="const RooAbsData&") {
	RooFactoryWSTool::as_DATA(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_DATA(%d)",i) ;
      } else if ((*ti)=="RooDataSet&" || (*ti)=="const RooDataSet&") {
	RooFactoryWSTool::as_DSET(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_DSET(%d)",i) ;
      } else if ((*ti)=="RooDataHist&" || (*ti)=="const RooDataHist&") {
	RooFactoryWSTool::as_DHIST(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_DHIST(%d)",i) ;
      } else if ((*ti)=="const RooArgSet&") {
	RooFactoryWSTool::as_SET(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_SET(%d)",i) ;
      } else if ((*ti)=="const RooArgList&") {
	RooFactoryWSTool::as_LIST(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_LIST(%d)",i) ;
      } else if ((*ti)=="const char*") {
	RooFactoryWSTool::as_STRING(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_STRING(%d)",i) ;	
      } else if ((*ti)=="Int_t" || (*ti)=="int" || (*ti)=="Bool_t" || (*ti)=="bool") {
	RooFactoryWSTool::as_INT(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_INT(%d)",i) ;	
      } else if ((*ti)=="Double_t") {
	RooFactoryWSTool::as_DOUBLE(i) ;
	cintExpr += Form(",RooFactoryWSTool::as_DOUBLE(%d)",i) ;	
      } else if (isEnum(ti->c_str())) {

	string qualvalue ;
	if (_args[i].find(Form("%s::",className)) != string::npos) {		    
	  qualvalue = _args[i].c_str() ;
	} else {	
	  qualvalue =  Form("%s::%s",className,_args[i].c_str()) ;	    
	}
	if (isValidEnumValue(ti->c_str(),qualvalue.c_str())) {
	  cintExpr += Form(",(%s)%s",ti->c_str(),qualvalue.c_str()) ;
	} else {
	  throw string(Form("Supplied argument %s does not represent a valid state of enum %s",_args[i].c_str(),ti->c_str())) ;
	  }
      } else {
	// Check if generic object store has argument of given name and type
	TObject& obj = RooFactoryWSTool::as_OBJ(i) ;

	// Strip argument type to bare type (i.e. const X& -> X)
	string btype ;
	if (ti->find("const ")==0) {
	  btype = ti->c_str()+6 ;
	} else {
	  btype = *ti ;
	}
	if (btype.find("&")) {
	  btype.erase(btype.size()-1,btype.size()) ;
	}

	// If btype if a typedef, substitute it by the true type name
	btype = string(TEnum::GetEnum(btype.c_str())->GetName());

	if (obj.InheritsFrom(btype.c_str())) {
	  cintExpr += Form(",(%s&)RooFactoryWSTool::as_OBJ(%d)",ti->c_str(),i) ;
	} else {
	  throw string(Form("Required argument with name %s of type '%s' is not in the workspace",_args[i].c_str(),ti->c_str())) ;
	}
      }
    }
    cintExpr += ") ;" ;
  } catch (const string &err) {
    coutE(ObjectHandling) << "RooFactoryWSTool::createArg() ERROR constructing " << className << "::" << objName << ": " << err << endl ;
    logError() ;
    return 0 ;
  }
  
  cxcoutD(ObjectHandling) << "RooFactoryWSTool::createArg() Construct expression is " << cintExpr << endl ;

  // Call CINT to perform constructor call. Catch any error thrown by argument conversion method
  RooAbsArg* arg = (RooAbsArg*) gROOT->ProcessLineFast(cintExpr.c_str()) ;

  if (arg) {    
    if (string(className)=="RooGenericPdf") {
      arg->setStringAttribute("factory_tag",Form("EXPR::%s(%s)",objName,varList)) ;
    } else if (string(className)=="RooFormulaVar") {
      arg->setStringAttribute("factory_tag",Form("expr::%s(%s)",objName,varList)) ;
    } else {
      arg->setStringAttribute("factory_tag",Form("%s::%s(%s)",className,objName,varList)) ;
    }
    if (_ws->import(*arg,Silence())) logError() ;
    RooAbsArg* ret = _ws->arg(objName) ;
    delete arg ;
    return ret ;
  } else {
    coutE(ObjectHandling) << "RooFactoryWSTool::createArg() ERROR in CINT constructor call to create object" << endl ;
    logError() ;
    return 0 ;
  }
}

////////////////////////////////////////////////////////////////////////////////

RooAddPdf* RooFactoryWSTool::add(const char *objName, const char* specList, Bool_t recursiveCoefs)
{
  // Spec list is of form a*A,b*B,c*C,D [ *d]

  RooArgList pdfList ;
  RooArgList coefList ;
  RooArgList pdfList2 ;

  try {

    char buf[BUFFER_SIZE] ;
    strlcpy(buf,specList,BUFFER_SIZE) ;
    char* save ;
    char* tok = R__STRTOK_R(buf,",",&save) ;
    while(tok) {
      char* star=strchr(tok,'*') ;
      if (star) {
	*star=0 ;
	pdfList.add(asPDF(star+1)) ;
	coefList.add(asFUNC(tok)) ;
      } else {
	pdfList2.add(asPDF(tok)) ;
      }
      tok = R__STRTOK_R(0,",",&save) ;
    }
    pdfList.add(pdfList2) ;

  } catch (const string &err) {
    coutE(ObjectHandling) << "RooFactoryWSTool::add(" << objName << ") ERROR creating RooAddPdf: " << err << endl ;    
    logError() ;
    return 0 ;
  }
  
  RooAddPdf* pdf =  new RooAddPdf(objName,objName,pdfList,coefList,recursiveCoefs) ;
  pdf->setStringAttribute("factory_tag",Form("SUM::%s(%s)",objName,specList)) ;
  if (_ws->import(*pdf,Silence())) logError() ;
  return (RooAddPdf*) _ws->pdf(objName) ;
}


////////////////////////////////////////////////////////////////////////////////

RooRealSumPdf* RooFactoryWSTool::amplAdd(const char *objName, const char* specList)
{
  // Spec list is of form a*A,b*B,c*C,D [ *d]

  RooArgList amplList ;
  RooArgList coefList ;
  RooArgList amplList2 ;

  try {

    char buf[BUFFER_SIZE] ;
    strlcpy(buf,specList,BUFFER_SIZE) ;
    char* save ;
    char* tok = R__STRTOK_R(buf,",",&save) ;
    while(tok) {
      char* star=strchr(tok,'*') ;
      if (star) {
	*star=0 ;
	amplList.add(asFUNC(star+1)) ;
	coefList.add(asFUNC(tok)) ;
      } else {
	amplList2.add(asFUNC(tok)) ;
      }
      tok = R__STRTOK_R(0,",",&save) ;
    }
    amplList.add(amplList2) ;

  } catch (const string &err) {
    coutE(ObjectHandling) << "RooFactoryWSTool::add(" << objName << ") ERROR creating RooRealSumPdf: " << err << endl ;    
    logError() ;
    return 0 ;
  }
  
  RooRealSumPdf* pdf =  new RooRealSumPdf(objName,objName,amplList,coefList,(amplList.getSize()==coefList.getSize())) ;
  pdf->setStringAttribute("factory_tag",Form("ASUM::%s(%s)",objName,specList)) ;
  if (_ws->import(*pdf,Silence())) logError() ;
  return (RooRealSumPdf*) _ws->pdf(objName) ;
}


////////////////////////////////////////////////////////////////////////////////

RooProdPdf* RooFactoryWSTool::prod(const char *objName, const char* pdfList) 
{
  _of = this ;

  // Separate conditional and non-conditional p.d.f terms
  RooLinkedList cmdList ;
  string regPdfList="{" ;
  char buf[BUFFER_SIZE] ;
  strlcpy(buf,pdfList,BUFFER_SIZE) ;
  char* save ;
  char* tok = R__STRTOK_R(buf,",",&save) ;
  while(tok) {
    char *sep = strchr(tok,'|') ;
    if (sep) {
      // Conditional term
      *sep=0 ;
      sep++ ;

      // |x is conditional on x, |~x is conditional on all but x
      Bool_t invCond(kFALSE) ;
      if (*sep=='~') {
	invCond=kTRUE ;
	sep++ ;
      }
      
      try {
 	cmdList.Add(Conditional(asSET(tok),asSET(sep),!invCond).Clone()) ;
      } catch (const string &err) {
	coutE(ObjectHandling) << "RooFactoryWSTool::prod(" << objName << ") ERROR creating RooProdPdf Conditional argument: " << err << endl ;
	logError() ;
	return 0 ;
      }
      
    } else {
      // Regular term
      if (regPdfList.size()>1) {
	regPdfList += "," ;
      }
      regPdfList += tok ;
    }
    tok = R__STRTOK_R(0,",",&save) ;
  }
  regPdfList += "}" ;
  
  RooProdPdf* pdf = 0 ;
  try {
    pdf = new RooProdPdf(objName,objName,asSET(regPdfList.c_str()),cmdList) ;
  } catch (const string &err) {
    coutE(ObjectHandling) << "RooFactoryWSTool::prod(" << objName << ") ERROR creating RooProdPdf input set of regular p.d.f.s: " << err << endl ;
    logError() ;
    pdf = 0 ;
  }
  cmdList.Delete() ;
  
  if (pdf) {
    pdf->setStringAttribute("factory_tag",Form("PROD::%s(%s)",objName,pdfList)) ;
    if (_ws->import(*pdf,Silence())) logError() ;
    delete pdf ;
    return (RooProdPdf*) _ws->pdf(objName) ;
  } else {
    return 0 ;
  }
}



////////////////////////////////////////////////////////////////////////////////

RooSimultaneous* RooFactoryWSTool::simul(const char* objName, const char* indexCat, const char* pdfMap)
{
  map<string,RooAbsPdf*> theMap ;
  // Add p.d.f. to index state mappings
  char buf[BUFFER_SIZE] ;
  strlcpy(buf,pdfMap,BUFFER_SIZE) ;
  char* save ;
  char* tok = R__STRTOK_R(buf,",",&save) ;
  while(tok) {
    char* eq = strchr(tok,'=') ;
    if (!eq) {
      coutE(ObjectHandling) << "RooFactoryWSTool::simul(" << objName << ") ERROR creating RooSimultaneous::" << objName 
			    << " expect mapping token of form 'state=pdfName', but found '" << tok << "'" << endl ;
      logError() ;
      return 0 ;
    } else {
      *eq = 0 ;

      try {
	theMap[tok] = &asPDF(eq+1) ;
      } catch (const string &err ) {
	coutE(ObjectHandling) << "RooFactoryWSTool::simul(" << objName << ") ERROR creating RooSimultaneous: " << err << endl ;
	logError() ;
      }
    }
    tok = R__STRTOK_R(0,",",&save) ;
  }


  // Create simultaneous p.d.f.
  RooSimultaneous* pdf(0) ;
  try {
    pdf = new RooSimultaneous(objName,objName,theMap,asCATLV(indexCat)) ;
  } catch (const string &err) {
    coutE(ObjectHandling) << "RooFactoryWSTool::simul(" << objName << ") ERROR creating RooSimultaneous::" << objName << " " << err << endl ;
    logError() ;
  }

  // Import p.d.f into workspace
  pdf->setStringAttribute("factory_tag",Form("SIMUL::%s(%s,%s)",objName,indexCat,pdfMap)) ;
  if (_ws->import(*pdf,Silence())) logError() ;
  return (RooSimultaneous*) _ws->pdf(objName) ;
}




////////////////////////////////////////////////////////////////////////////////

RooAddition* RooFactoryWSTool::addfunc(const char *objName, const char* specList) 
{
  RooArgList sumlist1 ;
  RooArgList sumlist2 ;

  try {

    char buf[BUFFER_SIZE] ;
    strlcpy(buf,specList,BUFFER_SIZE) ;
    char* save ;
    char* tok = R__STRTOK_R(buf,",",&save) ;
    while(tok) {
      char* star=strchr(tok,'*') ;
      if (star) {
	*star=0 ;
	sumlist2.add(asFUNC(star+1)) ;
	sumlist1.add(asFUNC(tok)) ;
      } else {
	sumlist1.add(asFUNC(tok)) ;
      }
      tok = R__STRTOK_R(0,",",&save) ;
    }

  } catch (const string &err) {
    coutE(ObjectHandling) << "RooFactoryWSTool::addfunc(" << objName << ") ERROR creating RooAddition: " << err << endl ;
    logError() ;
    return 0 ;
  }

  if (sumlist2.getSize()>0 && (sumlist1.getSize()!=sumlist2.getSize())) {
    coutE(ObjectHandling) << "RooFactoryWSTool::addfunc(" << objName << ") ERROR creating RooAddition: syntax error: either all sum terms must be products or none" << endl ;
    logError() ;
    return 0 ;
  }


  RooAddition* sum ;
  if (sumlist2.getSize()>0) {
    sum = new RooAddition(objName,objName,sumlist1,sumlist2) ;
  } else {
    sum = new RooAddition(objName,objName,sumlist1) ;
  }

  sum->setStringAttribute("factory_tag",Form("sum::%s(%s)",objName,specList)) ;
  if (_ws->import(*sum,Silence())) logError() ;
  delete sum ;
  return (RooAddition*) _ws->pdf(objName) ;
  
}




////////////////////////////////////////////////////////////////////////////////

RooProduct* RooFactoryWSTool::prodfunc(const char *objName, const char* pdfList) 
{
  return (RooProduct*) createArg("RooProduct",objName,Form("{%s}",pdfList)) ;
}





////////////////////////////////////////////////////////////////////////////////
/// Create a RooFit object from the given expression.
///
/// <table>
/// <tr><th> Creating variables <th>
/// <tr><td> `x[-10,10]`             <td>  Create variable x with given range and put it in workspace
/// <tr><td> `x[3,-10,10]`           <td>  Create variable x with given range and initial value and put it in workspace
/// <tr><td> `x[3]`                  <td>  Create variable x with given constant value
/// <tr><td> `<numeric literal>`     <td> Numeric literal expressions (0.5, -3 etc..) are converted to a RooConst(<numeric literal>)
///                                       wherever a RooAbsReal or RooAbsArg argument is expected
/// <tr><th> Creating categories <th>
/// <tr><td> `c[lep,kao,nt1,nt2]`    <td>  Create category c with given state names
/// <tr><td> `tag[B0=1,B0bar=-1]`    <td>  Create category tag with given state names and index assignments
/// <tr><th> Creating functions and p.d.f.s <th>
/// <tr><td> `MyPdf::g(x,m,s)`       <td> Create p.d.f or function of type MyPdf with name g with argument x,m,s
///                         Interpretation and number of arguments are mapped to the constructor arguments of the class 
///                         (after the name and title).
/// <tr><td> `MyPdf(x,m,s)`          <td> As above, but with an implicitly defined (unique) object name
/// <tr><th> Creating sets and lists (to be used as inputs above) <th>
/// <tr><td> `{a,b,c}`               <td> Create RooArgSet or RooArgList (as determined by context) from given contents
/// </table>
///
///
/// Objects that are not created, are assumed to exist in the workspace
/// Object creation expressions as shown above can be nested, e.g. one can do
/// ```
///   RooGaussian::g(x[-10,10],m[0],3) 
/// ```
/// to create a p.d.f and its variables in one go. This nesting can be applied recursively e.g.
/// ```
///   SUM::model( f[0.5,0,1] * RooGaussian::g( x[-10,10], m[0], 3] ),
///                            RooChebychev::c( x, {a0[0.1],a1[0.2],a2[-0.3]} ))
/// ```
/// creates the sum of a Gaussian and a Chebychev and all its variables.
///
///
/// A seperate series of operator meta-type exists to simplify the construction of composite expressions
/// meta-types in all capitals (SUM) create p.d.f.s, meta types in lower case (sum) create
/// functions.
///
/// <table>
/// <tr><th> Expression <th> Effect
/// <tr><td> `SUM::name(f1*pdf1,f2*pdf2,pdf3]`  <td> Create sum p.d.f name with value f1*pdf1+f2*pdf2+(1-f1-f2)*pdf3
/// <tr><td> `RSUM::name(f1*pdf1,f2*pdf2,pdf3]` <td> Create recursive sum p.d.f. name with value f1*pdf1 + (1-f1)(f2*pdf2 + (1-f2)pdf3)
/// <tr><td> `ASUM::name(f1*amp1,f2*amp2,amp3]` <td> Create sum p.d.f. name with value f1*amp1+f2*amp2+(1-f1-f2)*amp3 where amplX are amplitudes of type RooAbsReal
/// <tr><td> `sum::name(a1,a2,a3]`              <td> Create sum function with value a1+a2+a3
/// <tr><td> `sum::name(a1*b1,a2*b2,a3*b 3]`    <td> Create sum function with value a1*b1+a2*b2+a3*b3
/// <tr><td> `PROD::name(pdf1,pdf2]`            <td> Create product of p.d.f with 'name' with given input p.d.fs
/// <tr><td> `PROD::name(pdf1|x,pdf2]`          <td> Create product of conditional p.d.f. pdf1 given x and pdf2
/// <tr><td> `prod::name(a,b,c]`                <td> Create production function with value a*b*c
/// <tr><td> `SIMUL::name(cat,a=pdf1,b=pdf2]`   <td> Create simultaneous p.d.f index category cat. Make pdf1 to state a, pdf2 to state b
/// <tr><td> `EXPR::name(<expr>,var,...]`       <td> Create a generic p.d.f that interprets the given expression
/// <tr><td> `expr::name(<expr>,var,...] `       <td> Create a generic function that interprets the given expression
/// </table>
///
/// The functionality of high-level object creation tools like RooSimWSTool, RooCustomizer and RooClassFactory
/// is also interfaced through meta-types in the factory.
/// <table>
/// <tr><th> Interface to %RooSimWSTool <th>
/// <tr><td> `SIMCLONE::name( modelPdf, $ParamSplit(...), $ParamSplitConstrained(...), $Restrict(...) ]`
///             <td> Clone-and-customize modelPdf according to ParamSplit and ParamSplitConstrained()
///                  specifications and return a RooSimultaneous p.d.f. of all built clones
///
/// <tr><td> `MSIMCLONE::name( masterIndex, $AddPdf(mstate1, modelPdf1, $ParamSplit(...)), $AddPdf(mstate2,modelPdf2),...) ]`
///                        <td> Clone-and-customize multiple models (modelPdf1,modelPdf2) according to ParamSplit and
///                                                                             ParamSplitConstrained() specifications and return a RooSimultaneous p.d.f. of all built clones,
///                                                                             using the specified master index to map prototype p.d.f.s to master states
/// <tr><th> Interface to %RooCustomizer <th>
/// <tr><td> `EDIT::name( orig, substNode=origNode), ... ]`                             <td> Create a clone of input object orig, with the specified replacements operations executed
/// <tr><td> `EDIT::name( orig, origNode=$REMOVE(), ... ]`                              <td> Create clone of input removing term origNode from all PROD() terms that contained it
/// <tr><td> `EDIT::name( orig, origNode=$REMOVE(prodname,...), ... ]`                  <td> As above, but restrict removal of origNode to PROD term(s) prodname,...
///
///
/// <tr><th> Interface to %RooClassFactory <th>
/// <tr><td> `CEXPR::name(<expr>,var,...]`       <td> Create a custom compiled p.d.f that evaluates the given expression
/// <tr><td> `cexpr::name(<expr>,var,...]`       <td> Create a custom compiled function that evaluates the given expression
///
///
/// <tr><td> `$MetaType(...)`        <td> Meta argument that does not result in construction of an object but is used logically organize
///                         input arguments in certain operator p.d.f. constructions. The defined meta arguments are context dependent.
///                         The only meta argument that is defined globally is `$Alias(typeName,aliasName)` to
///                         define aliases for type names. For the definition of meta arguments in operator p.d.f.s
///                         see the definitions below.
/// </table>
RooAbsArg* RooFactoryWSTool::process(const char* expr) 
{

//   cout << "RooFactoryWSTool::process() " << expr << endl ;

  // First perform basic syntax check
  if (checkSyntax(expr)) {
    return 0 ;
  }

  // Allocate work buffer
  char* buf = new char[strlen(expr)+1] ;

  // Copy to buffer while absorbing white space and newlines
  char* buftmp = buf ;
  while(*expr) {
    if (!isspace(*expr)) {
      *buftmp = *expr ;
      buftmp++ ;
    }
    expr++ ;
  }
  *buftmp=0 ;


  // Clear error count and start a transaction in the workspace
  clearError() ;
  ws().startTransaction() ;

  // Process buffer
  string out ;
  try {
    out = processExpression(buf) ;
  } catch (const string &error) {
    coutE(ObjectHandling) << "RooFactoryWSTool::processExpression() ERROR in parsing: " << error << endl ;
    logError() ;
  }

  // If there were no errors commit the transaction, cancel it otherwise
  if (errorCount()>0) {
    coutE(ObjectHandling) << "RooFactoryWSTool::processExpression() ERRORS detected, transaction to workspace aborted, no objects committed" << endl ;
    ws().cancelTransaction() ;
  } else {
    ws().commitTransaction() ;
  }

  
  // Delete buffer 
  delete[] buf ;

  return out.size() ? ws().arg(out.c_str()) : 0 ;
}




////////////////////////////////////////////////////////////////////////////////
/// Process a single high-level expression or list of
/// expressions. The returned string a the reduced expression where
/// all inline object creations have been executed and substituted
/// with the name of the created object
///
/// e.g. 'RooGaussian::g(x,m,s)' --> 'g'
///      '{x(-10,10),s}          --> '{x,s}'

std::string RooFactoryWSTool::processExpression(const char* token) 
{
  // Delegate handling to list processor if token starts with {, otherwise
  // call single expression processor
  if (string(token).find("$Alias(")==0) {
    processAliasExpression(token) ;
  }

  if (token[0]=='{') {    
    // Process token as list if it starts with '{'
    return processListExpression(token) ;
  } else {
    // Process token as single item otherwise
    return processCompositeExpression(token) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Process a single composite expression 
///
/// e.g. 'A=RooGaussian::g[x,m,s]' --> 'A=g'
/// e.g. 'f[0,1]*RooGaussian::g[x,m,s]' --> 'f*g'
/// e.g. 'RooGaussian::g(x,y,s)|x' --> g|x'
/// e.g. '$MetaArg(RooGaussian::g[x,m,s],blah)' --> '$MetaArg(g,blah)'

std::string RooFactoryWSTool::processCompositeExpression(const char* token) 
{
  // Allocate and fill work buffer
   const size_t bufBaseSize = strlen(token)+1;
  char* buf_base = new char[bufBaseSize] ;
  char* buf = buf_base ;
  strlcpy(buf,token,bufBaseSize) ;
  char* p = buf ;

  list<string> singleExpr ;
  list<char> separator ;
  Int_t blevel(0) ;
  Bool_t litmode(kFALSE) ;
  while(*p) {

    // Keep track of opening and closing brackets
    if (*p=='{' || *p=='(' || *p=='[') blevel++ ;
    if (*p=='}' || *p==')' || *p==']') blevel-- ;
    
    // Keep track of string literals
    if (*p=='"' || *p=='\'') litmode = !litmode ;

    // If we are zero-bracket level and encounter a |, store
    // the remainder of the string as suffix and exit loop
    if (!litmode && blevel==0 && ( (*p)=='=' || (*p) == '|' || (*p) == '*')) {
      separator.push_back(*p) ;
      *p=0 ;
      singleExpr.push_back(buf) ;
      buf = p+1 ;
    }
    p++ ;
  }
  if (*buf) {
    singleExpr.push_back(buf) ;
  }
  if (singleExpr.size()==1) {    
    string ret = processSingleExpression(token) ;
    delete[] buf_base ;
    return ret ;
  }

  string ret ;
  list<char>::iterator ic = separator.begin() ;
  for (list<string>::iterator ii = singleExpr.begin() ; ii!=singleExpr.end() ; ++ii) {
    ret += processSingleExpression(ii->c_str()) ;
    if (ic != separator.end()) {
      ret += *ic ;
      ++ic ;
    }
  }

  delete[] buf_base ;
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Process a single high-level expression. The returned string a the reduced
/// expression where all inline object creations have been executed and substituted
/// with the name of the created object
///
/// e.g. 'RooGaussian::g(x,m,s)' --> 'g'
/// e.g. 'x[-10,10]' --> 'x'

std::string RooFactoryWSTool::processSingleExpression(const char* arg) 
{
  // Handle empty strings here
  if (strlen(arg)==0) {
    return string("") ;
  }

  // Handle string literal case
  if (arg[0]=='\'' || arg[0]=='"') {
    return string(arg) ;
  }

  // Allocate and fill work buffer
  const size_t bufSize = strlen(arg)+1;
  char* buf = new char[bufSize] ;
  strlcpy(buf,arg,bufSize) ;
  char* bufptr = buf ;

  string func,prefix ;
  vector<string> args ;

  // Process token into arguments
  char* save ;
  char* tmpx = R__STRTOK_R(buf,"([",&save) ;
  func = tmpx ? tmpx : "" ;
  char* p = R__STRTOK_R(0,"",&save) ;
  
  // Return here if token is fundamental
  if (!p) {
    delete[] buf ;
    return arg ;
  }
    

  char* tok = p ;
  Int_t blevel=0 ;
  Bool_t litmode(kFALSE) ;
  while(*p) {

    // Keep track of opening and closing brackets
    if (*p=='{' || *p=='(' || *p=='[') blevel++ ;
    if (*p=='}' || *p==')' || *p==']') blevel-- ;

    // Keep track of string literals
    if (*p=='"' || *p=='\'') litmode = !litmode ;


    // If we encounter a comma at zero bracket level
    // finalize the current token as a completed argument
    // and start the next token
    if (!litmode && blevel==0 && ((*p)==',')) {
      *p = 0 ;
      args.push_back(tok) ;
      tok = p+1 ;
    }

    p++ ;
  }

  // If the last character was a closing bracket, kill
  // it in the buffer
  if (p>bufptr && (*(p-1)==')'||*(p-1)==']')) {
    *(p-1)=0 ;
  }

  // Finalize last token as argument
  string tmp = tok ;
  
  // If there is a suffix left in the work buffer attach it to 
  // this argument
  p = R__STRTOK_R(0,"",&save) ;
  if (p) tmp += p ;
  args.push_back(tmp) ;

  // Delete the work buffer
  delete[] buf ;

  // If function contains :: then call createArg to process this arg, otherwise
  // call createVariable
  string ret ;

  // Determine type of leading bracket
  char lb = ' ' ;  
  for(const char* pp=arg ; *pp!=0 ; pp++) {
    if (*pp=='(' || *pp=='[' || *pp=='{') {
      lb = *pp ;
      break ;
    }
  }

  if (strstr(func.c_str(),"::")) {
    if (lb=='(') {
      // Create function argument with instance name
      ret= processCreateArg(func,args) ;
    } else {
      coutE(ObjectHandling) << "RooFactoryWSTool::processSingleExpression(" << arg << "): ERROR: Syntax error: Class::Instance must be followed by (...)" << endl ;
      logError() ;
    }
  } else if (func[0]!='$'){    
    if (lb=='[') {
      // Create variable argument
      ret= processCreateVar(func,args) ;
    } else if (lb=='(') {

      // Create function argument with autoname
      string autoname ;
      if (!_autoNamePrefix.empty()) {
	// If we're inside a function creation call to a higher level object, use its
	// name as base for the autoname
	autoname = (Form("%s::%s",func.c_str(),_autoNamePrefix.top().c_str())) ;
      } else {
	// Otherwise find a free global_%d name
	static Int_t globCounter = 0 ;
	while(true) {
	  autoname = Form("gobj%d",globCounter) ;
	  globCounter++ ;
	  if (!ws().arg(autoname.c_str())) {
	    break ;
	  }
	}
	autoname = Form("%s::%s",func.c_str(),autoname.c_str()) ;
      }
      ret= processCreateArg(autoname,args) ;
    } else {
      coutE(ObjectHandling) << "RooFactoryWSTool::processSingleExpression(" << arg << "): ERROR: Syntax error: expect either Class(...) or Instance[...]" << endl ;
      logError() ;
    }
  } else {
    if (lb=='(') {
      // Process meta function (compile arguments, but not meta-function itself)
      ret= processMetaArg(func,args) ;
    } else {
      coutE(ObjectHandling) << "RooFactoryWSTool::processSingleExpression(" << arg << "): ERROR: Syntax error: $MetaClass must be followed by (...)" << endl ;
      logError() ;
    }
  }

  // Return reduced token with suffix
  return ret ;
}


////////////////////////////////////////////////////////////////////////////////
/// Process a list of high-level expression. The returned string a the reduced
/// expression list where all inline object creations have been executed and substituted
/// with the name of the created object
///
/// E.g.   '{x(-10,10),s}  --> '{x,s}'

string RooFactoryWSTool::processListExpression(const char* arg) 
{
  // Allocate and fill work buffer
  const size_t bufSize = strlen(arg)+1;
  char* buf = new char[bufSize] ;
  strlcpy(buf,arg,bufSize) ;

  vector<string> args ;

  // Start running pointer at position 1 to skip opening bracket
  char* tok = buf+1 ;
  char* p = buf+1 ;
  
  // Processing look
  Int_t level(0) ;
  while(*p) {

    // Track bracketing level
    if (*p=='{' || *p=='(' || *p=='[') level++ ;
    if (*p=='}' || *p==')' || *p==']') level-- ;

    
    // If we encounter a comma at zero bracket level
    // finalize the current token as a completed argument
    // and start the next token
    if (level==0 && ((*p)==',')) {
      *p = 0 ;
      args.push_back(tok) ;
      tok = p+1 ;
    }

    p++ ;
  }

  // Finalize token as last argument
  if (p>buf && *(p-1)=='}') {
    *(p-1)=0 ;
  }
  args.push_back(tok) ;

  // Delete work buffer
  delete[] buf ;

  // Process each argument in list and construct reduced 
  // expression to be returned
  string ret("{") ;
  vector<string>::iterator iter = args.begin() ;
  Int_t i(0) ;
  while(iter!= args.end()) {
    if (strlen(ret.c_str())>1) ret += "," ;
    if (!_autoNamePrefix.empty()) {
      _autoNamePrefix.push(Form("%s%d",_autoNamePrefix.top().c_str(),i+1)) ;
    }
    ret += processSingleExpression(iter->c_str()) ;
    if (!_autoNamePrefix.empty()) {
      _autoNamePrefix.pop() ;
    }
    ++iter ;
    i++ ;
  }
  ret += "}" ;

  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Parse token

string RooFactoryWSTool::processAliasExpression(const char* token) 
{
  vector<string> args = splitFunctionArgs(token) ;
  if (args.size()!=2) {
    coutE(ObjectHandling) << "RooFactorWSTool::processAliasExpression() ERROR $Alias() takes exactly two arguments, " << args.size() << " args found" << endl ;
    logError() ;
    return string() ;
  }

  // Insert alias in table
  _typeAliases[args[1]] = args[0] ;
  
  return string() ;
}




////////////////////////////////////////////////////////////////////////////////

TClass* RooFactoryWSTool::resolveClassName(const char* className) 
{
  // First do recursive alias expansion
  while (true) {
    map<string,string>::iterator item = _typeAliases.find(className) ;
    
    // If an alias is found, recurse
    if (item != _typeAliases.end()) {
      className = item->second.c_str() ;
    } else {
      break ;
    }
  }

  // Now find dealiased class in ROOT class table
  TClass* tc =  TClass::GetClass(className,kTRUE,kTRUE) ;

  // If its not there, try prefixing with Roo
  if (!tc) {
    tc = TClass::GetClass(Form("Roo%s",className)) ;
    if (!tc) {
      coutE(ObjectHandling) << "RooFactoryWSTool::createArg() ERROR class " << className << " not defined in ROOT class table" << endl ;
      logError() ;
      return 0 ;
    } 
  }
  return tc ;
}



////////////////////////////////////////////////////////////////////////////////

string RooFactoryWSTool::varTag(string& func, vector<string>& args)
{
  string ret ;
  ret += func ;
  ret += "[" ;
  for (vector<string>::iterator iter = args.begin() ; iter!=args.end() ; ++iter) {
    if (iter!=args.begin()) {
      ret += "," ;
    }
    ret += *iter ;
  }
  ret += "]" ;
  return ret ;
}




////////////////////////////////////////////////////////////////////////////////
/// Glue function between high-level syntax and low-level factory call to createVariable:
/// Process a parsed call to create a variable named 'func'
/// 
/// If initial token is non-numeric, a RooCategory will be created, and the args are interpreted
/// as either state names or 'name=id' assignments. Otherwise a RooRealvar is created and the
/// arg list is interpreted as follows:
/// If list has two args, these are interpreted as xmin,xmax
/// If list has three args, these are interpreted as xinit,xmin,xmax
/// If list has one arg, this is interpreted as xinit and the variable is set as constant

string RooFactoryWSTool::processCreateVar(string& func, vector<string>& args)
{  

  // Determine if first arg is numeric
  string first = *(args.begin()) ;
  if (isdigit(first[0]) || first[0]=='.' || first[0]=='+' || first[0]=='-') {
    
    // Create a RooRealVar
    vector<string>::iterator ai = args.begin() ;
    if (args.size()==1) {

      // One argument, create constant variable with given value
      Double_t xinit = atof((ai)->c_str()) ;
      cxcoutD(ObjectHandling) << "CREATE variable " << func << " xinit = " << xinit << endl ;
      RooRealVar tmp(func.c_str(),func.c_str(),xinit) ;
      tmp.setStringAttribute("factory_tag",varTag(func,args).c_str()) ;
      if (_ws->import(tmp,Silence())) {
	logError() ;
      }

    } else if (args.size()==2) {

      // Two arguments, create variable with given range
      Double_t xlo = atof((ai++)->c_str()) ;
      Double_t xhi = atof(ai->c_str()) ;
      cxcoutD(ObjectHandling) << "CREATE variable " << func << " xlo = " << xlo << " xhi = " << xhi << endl ;
      RooRealVar tmp(func.c_str(),func.c_str(),xlo,xhi) ;
      tmp.setStringAttribute("factory_tag",varTag(func,args).c_str()) ;
      if (_ws->import(tmp,Silence())) {
	logError() ;
      }

    } else if (args.size()==3) {

      // Three arguments, create variable with given initial value and range
      Double_t xinit = atof((ai++)->c_str()) ;
      Double_t xlo = atof((ai++)->c_str()) ;
      Double_t xhi = atof(ai->c_str()) ;
      cxcoutD(ObjectHandling) << "CREATE variable " << func << " xinit = " << xinit << " xlo = " << xlo << " xhi = " << xhi << endl ;
      RooRealVar tmp(func.c_str(),func.c_str(),xinit,xlo,xhi) ;
      tmp.setStringAttribute("factory_tag",varTag(func,args).c_str()) ;
      if (_ws->import(tmp,Silence())) {
	logError() ;
      }
    }
  } else {

    // Create a RooAbsCategory
    string allStates ;
    for (vector<string>::iterator ai = args.begin() ; ai!=args.end() ; ++ai) {
      if (allStates.size()>0) {
	allStates += "," ;
      }
      allStates += *ai ;
    }
    createCategory(func.c_str(),allStates.c_str()) ;

  }
  return func ;
}
      

////////////////////////////////////////////////////////////////////////////////
/// Glue function between high-level syntax and low-level factory call to createArg:
/// Process a parsed call to create a p.d.f named func
///
/// The func arg is interpreted as ClassName::ObjectName and the arglist is passed
/// verbatim to createArg. The received arglist is expected to be fully reduced (i.e.
/// all inline object creations must have been compiled)

string RooFactoryWSTool::processCreateArg(string& func, vector<string>& args)
{  
  // Allocate and fill work buffer
  char buf[BUFFER_SIZE] ;
  strlcpy(buf,func.c_str(),BUFFER_SIZE) ;

  // Split function part in class name and instance name
  char* save ;
  const char *className = R__STRTOK_R(buf,":",&save) ;
  const char *instName = R__STRTOK_R(0,":",&save) ;
  if (!className) className = "";
  if (!instName) instName = "" ;

  // Concatenate list of args into comma separated string
  char pargs[BUFFER_SIZE] ;
  pargs[0] = 0 ;
  vector<string>::iterator iter = args.begin() ;
  vector<string> pargv ;
  Int_t iarg(0) ;
  while(iter!=args.end()) {
    if (strlen(pargs)>0) strlcat(pargs,",",BUFFER_SIZE) ;
    _autoNamePrefix.push(Form("%s_%d",instName,iarg+1)) ;
    string tmp = processExpression(iter->c_str()) ;
    _autoNamePrefix.pop() ;
    strlcat(pargs,tmp.c_str(),BUFFER_SIZE) ;
    pargv.push_back(tmp) ;
    ++iter ;
    iarg++ ;
  }

  // Look up if func is a special
  for (map<string,IFace*>::iterator ii=hooks().begin() ; ii!=hooks().end() ; ++ii) {
  }
  if (hooks().find(className) != hooks().end()) {
    IFace* iface = hooks()[className] ;
    return iface->create(*this, className,instName,pargv) ;
  }

  createArg(className,instName,pargs) ;

  return string(instName) ;
}
      


////////////////////////////////////////////////////////////////////////////////
/// Concatenate list of args into comma separated string

std::string RooFactoryWSTool::processMetaArg(std::string& func, std::vector<std::string>& args) 
{
  char pargs[BUFFER_SIZE] ;
  pargs[0] = 0 ;
  vector<string>::iterator iter = args.begin() ;
  vector<string> pargv ;
  while(iter!=args.end()) {
    if (strlen(pargs)>0) strlcat(pargs,",",BUFFER_SIZE) ;
    string tmp = processExpression(iter->c_str()) ;
    strlcat(pargs,tmp.c_str(),BUFFER_SIZE) ;
    pargv.push_back(tmp) ;
    ++iter ;
  }

  string ret = func+"("+pargs+")" ;  
  return ret ;
}




////////////////////////////////////////////////////////////////////////////////
/// Allocate and fill work buffer

vector<string> RooFactoryWSTool::splitFunctionArgs(const char* funcExpr) 
{
  const size_t bufSize = strlen(funcExpr)+1;
  char* buf = new char[bufSize] ;
  strlcpy(buf,funcExpr,bufSize) ;
  char* bufptr = buf ;

  string func ;
  vector<string> args ;

  // Process token into arguments
  char* save ;
  char* tmpx = R__STRTOK_R(buf,"(",&save) ;  
  func = tmpx ? tmpx : "" ;
  char* p = R__STRTOK_R(0,"",&save) ;
  
  // Return here if token is fundamental
  if (!p) {
    delete[] buf ;
    return args ;
  }

  char* tok = p ;
  Int_t blevel=0 ;
  Bool_t litmode(kFALSE) ;
  while(*p) {

    // Keep track of opening and closing brackets
    if (*p=='{' || *p=='(' || *p=='[') blevel++ ;
    if (*p=='}' || *p==')' || *p==']') blevel-- ;

    // Keep track of string literals
    if (*p=='"' || *p=='\'') litmode = !litmode ;


    // If we encounter a comma at zero bracket level
    // finalize the current token as a completed argument
    // and start the next token
    if (!litmode && blevel==0 && ((*p)==',')) {
      *p = 0 ;
      args.push_back(tok) ;
      tok = p+1 ;
    }

    p++ ;
  }

  // If the last character was a closing bracket, kill
  // it in the buffer
  if (p>bufptr && *(p-1)==')') {
    *(p-1)=0 ;
  }

  // Finalize last token as argument
  string tmp = tok ;
  
  // If there is a suffix left in the work buffer attach it to 
  // this argument
  p = R__STRTOK_R(0,"",&save) ;
  if (p) tmp += p ;
  args.push_back(tmp) ;

  // Delete the work buffer
  delete[] buf ;
  
  return args ;
}





////////////////////////////////////////////////////////////////////////////////
/// Perform basic syntax on given factory expression. If function returns
/// true syntax errors are found.

Bool_t RooFactoryWSTool::checkSyntax(const char* arg) 
{
  // Count parentheses
  Int_t nParentheses(0), nBracket(0), nAccolade(0) ;
  const char* ptr = arg ;
  while(*ptr) {
    if (*ptr=='(') nParentheses++ ;
    if (*ptr==')') nParentheses-- ;
    if (*ptr=='[') nBracket++ ;
    if (*ptr==']') nBracket-- ;
    if (*ptr=='{') nAccolade++ ;
    if (*ptr=='}') nAccolade-- ;
    ptr++ ;
  }
  if (nParentheses!=0) {
    coutE(ObjectHandling) << "RooFactoryWSTool::checkSyntax ERROR non-matching '" << (nParentheses>0?"(":")") << "' in expression" << endl ;
    return kTRUE ;
  }
  if (nBracket!=0) {
    coutE(ObjectHandling) << "RooFactoryWSTool::checkSyntax ERROR non-matching '" << (nBracket>0?"[":"]") << "' in expression" << endl ;
    return kTRUE ;
  }
  if (nAccolade!=0) {
    coutE(ObjectHandling) << "RooFactoryWSTool::checkSyntax ERROR non-matching '" << (nAccolade>0?"{":"}") << "' in expression" << endl ;
    return kTRUE ;
  }
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////

void RooFactoryWSTool::checkIndex(UInt_t idx) 
{
  if (idx>_of->_args.size()-1) {
    throw string(Form("Need argument number %d, but only %d args are provided",idx,(Int_t)_of->_args.size())) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooAbsArg reference found in workspace

RooAbsArg& RooFactoryWSTool::asARG(const char* arg) 
  {
  // If arg is a numeric string, make a RooConst() of it here
  if (arg[0]=='.' || arg[0]=='+' || arg[0] == '-' || isdigit(arg[0])) {
    return RooConst(atof(arg)) ;
  }

  // Otherwise look it up by name in the workspace
  RooAbsArg* rarg = ws().arg(arg) ;
  if (!rarg) {
    throw string(Form("RooAbsArg named %s not found",arg)) ;
  }
  return *rarg ;
}



////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooAbsReal reference found in workspace

RooAbsReal& RooFactoryWSTool::asFUNC(const char* arg) 
{
  // If arg is a numeric string, make a RooConst() of it here
  if (arg[0]=='.' || arg[0]=='+' || arg[0] == '-' || isdigit(arg[0])) {
    return RooConst(atof(arg)) ;
  }

  RooAbsArg* rarg = ws().arg(arg) ;
  if (!rarg) {
    throw string(Form("RooAbsReal named %s not found",arg)) ;
  }
  RooAbsReal* real = dynamic_cast<RooAbsReal*>(rarg) ;
  if (!real) {
    throw string(Form("Object named %s is not of type RooAbsReal",arg)) ;
  }
  return *real ;
}



////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooAbsRealLValue reference found in workspace

RooAbsRealLValue& RooFactoryWSTool::asVARLV(const char* arg) 
{
  // If arg is a numeric string, throw error as lvalue is required
  if (arg[0]=='.' || arg[0]=='+' || arg[0] == '-' || isdigit(arg[0])) {
    throw string(Form("Numeric literal provided for argument (%s), but lvalue is required",arg)) ;
  }

  RooAbsArg* rarg = ws().arg(arg) ;
  if (!rarg) {
    throw string(Form("RooAbsRealLValue named %s not found",arg)) ;
  }
  RooAbsRealLValue* reallv = dynamic_cast<RooAbsRealLValue*>(rarg) ;
  if (!reallv) {
    throw string(Form("Object named %s is not of type RooAbsRealLValue",arg)) ;
  }
  return *reallv ;
}



////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooRealVar reference found in workspace

RooRealVar& RooFactoryWSTool::asVAR(const char* arg) 
{
  RooRealVar* var = ws().var(arg) ;
  if (!var) {
    throw string(Form("RooRealVar named %s not found",arg)) ;
  }
  return *var ;
}




////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooAbsPdf reference found in workspace

RooAbsPdf& RooFactoryWSTool::asPDF(const char* arg) 
{
  RooAbsPdf* pdf = ws().pdf(arg) ;
  if (!pdf) {
    throw string(Form("RooAbsPdf named %s not found",arg)) ;
  }
  return *pdf ;
}




////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooResolutionModel reference found in workspace

RooResolutionModel& RooFactoryWSTool::asRMODEL(const char* arg) 
{
  RooAbsArg* rarg = ws().arg(arg) ;
  if (!rarg) {
    throw string(Form("RooResolutionModel named %s not found",arg)) ;
  }
  RooResolutionModel * rmodel = dynamic_cast<RooResolutionModel*>(rarg) ;
  if (!rmodel) {
    throw string(Form("Object named %s is not of type RooResolutionModel",arg)) ;
  }
  return *rmodel ;
}




////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooAbsCategory reference found in workspace

RooAbsCategory& RooFactoryWSTool::asCATFUNC(const char* arg) 
{
  RooAbsArg* rarg = ws().arg(arg) ;
  if (!rarg) {
    throw string(Form("RooAbsCategory named %s not found",arg)) ;
  }
  RooAbsCategory* catf = dynamic_cast<RooAbsCategory*>(rarg) ;
  if (!catf) {
    throw string(Form("Object named %s is not of type RooAbsCategory",arg)) ;
  }
  return *catf ;
}



////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooAbsCategoryLValue reference found in workspace

RooAbsCategoryLValue& RooFactoryWSTool::asCATLV(const char* arg) 
{
  RooAbsArg* rarg = ws().arg(arg) ;
  if (!rarg) {
    throw string(Form("RooAbsCategoryLValue named %s not found",arg)) ;
  }

  RooAbsCategoryLValue* catlv = dynamic_cast<RooAbsCategoryLValue*>(rarg) ;
  if (!catlv) {
    throw string(Form("Object named %s is not of type RooAbsCategoryLValue",arg)) ;
  }
  return *catlv ;
}



////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooCategory reference found in workspace

RooCategory& RooFactoryWSTool::asCAT(const char* arg) 
{
  RooCategory* cat = ws().cat(arg) ;
  if (!cat) {
    throw string(Form("RooCategory named %s not found",arg)) ;
  }
  return *cat ;
}





////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooArgSet of objects found in workspace

RooArgSet RooFactoryWSTool::asSET(const char* arg) 
{
  char tmp[BUFFER_SIZE] ;
  strlcpy(tmp,arg,BUFFER_SIZE) ;

  RooArgSet s ;
  
  // If given object is not of {,,,} form, interpret given string as name of defined set
  if (arg[0]!='{') {
    // cout << "asSet(arg='" << arg << "') parsing as defined set" << endl ;
    const RooArgSet* defSet = ws().set(arg) ;
    if (defSet) {
      // cout << "found defined set: " << *defSet << endl ;
      s.add(*defSet) ;
      return s ;
    }
  }

  char* save ;
  char* tok = R__STRTOK_R(tmp,",{}",&save) ;
  int i(0);
  while(tok) {

    // If arg is a numeric string, make a RooConst() of it here
    if (tok[0]=='.' || tok[0]=='+' || tok[0] == '-' || isdigit(tok[0])) {
      s.add(RooConst(atof(tok))) ;
    } else if (tok[0] == '\'') {
       tok[strlen(tok) - 1] = 0;
       RooStringVar *sv = new RooStringVar(Form("string_set_item%03d", i++), "string_set_item", tok + 1);
       s.add(*sv);
    } else {
      RooAbsArg* aarg = ws().arg(tok) ;
      if (aarg) {
	s.add(*aarg)  ;
      } else {
	throw string(Form("RooAbsArg named %s not found",tok)) ;
      }
    }
    tok = R__STRTOK_R(0,",{}",&save) ;
  }

  return s ;
}



////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooArgList of objects found in workspace

RooArgList RooFactoryWSTool::asLIST(const char* arg) 
{
  char tmp[BUFFER_SIZE] ;
  strlcpy(tmp,arg,BUFFER_SIZE) ;

  RooArgList l ;
  char* save ;
  char* tok = R__STRTOK_R(tmp,",{}",&save) ;
  while(tok) {

    // If arg is a numeric string, make a RooConst() of it here
    if (tok[0]=='.' || tok[0]=='+' || tok[0] == '-' || isdigit(tok[0])) {
      l.add(RooConst(atof(tok))) ;
    } else if (tok[0] == '\'') {
       tok[strlen(tok) - 1] = 0;
       RooStringVar *sv = new RooStringVar("listarg", "listarg", tok + 1);
       l.add(*sv);
    } else {
      RooAbsArg* aarg = ws().arg(tok) ;
      if (aarg) {
	l.add(*aarg)  ;
      } else {
	throw string(Form("RooAbsArg named %s not found",tok)) ;
      }
    }
    tok = R__STRTOK_R(0,",{}",&save) ;
  }
  
  return l ;
}



////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooAbsData object found in workspace

RooAbsData& RooFactoryWSTool::asDATA(const char* arg) 
{
  RooAbsData* data = ws().data(arg) ;
  if (!data) {
      throw string(Form("RooAbsData named %s not found",arg)) ;    
  }
  return *data ;
}



////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooDataHist object found in workspace

RooDataHist& RooFactoryWSTool::asDHIST(const char* arg) 
{
  RooAbsData* data = ws().data(arg) ;
  if (!data) {
    throw string(Form("RooAbsData named %s not found",arg)) ;    
  }
  RooDataHist* hist = dynamic_cast<RooDataHist*>(data) ;
  if (!hist) {
    throw string(Form("Dataset named %s is not of type RooDataHist",arg)) ;    
  }
  return *hist ;
}


////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as RooDataSet object found in workspace

RooDataSet& RooFactoryWSTool::asDSET(const char* arg) 
{
  RooAbsData* data = ws().data(arg) ;
  if (!data) {
    throw string(Form("RooAbsData named %s not found",arg)) ;    
  }
  RooDataSet* dset = dynamic_cast<RooDataSet*>(data) ;
  if (!dset) {
    throw string(Form("Dataset named %s is not of type RooDataSet",arg)) ;    
  }
  return *dset ;
}



////////////////////////////////////////////////////////////////////////////////

TObject& RooFactoryWSTool::asOBJ(const char* arg)
{
  TObject* obj = ws().obj(arg) ;
  if (!obj) {
    throw string(Form("Object named %s not found",arg)) ;    
  }
  return *obj ;
}



////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as const char* 

const char* RooFactoryWSTool::asSTRING(const char* arg) 
{
  static vector<string> cbuf(10) ;
  static unsigned int cbuf_idx = 0 ;

  // Handle empty string case: return null pointer
  if (arg==0 || strlen(arg)==0) {
    return 0 ;
  }
  
  // Fill cyclical buffer entry with quotation marked stripped version of string literal
  // and return pointer to stripped buffer
  cbuf[cbuf_idx].clear() ;
  const char* p = arg+1 ;
  while(*p && (*p) != '"' && (*p) !='\'' ) {
    cbuf[cbuf_idx] += *(p++) ;
  }  
  const char* ret = cbuf[cbuf_idx].c_str() ;

  // Increment buffer pointer by one
  cbuf_idx++ ;
  if (cbuf_idx==cbuf.size()) cbuf_idx=0 ;

  return ret ;
}


////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as Int_t

Int_t RooFactoryWSTool::asINT(const char* arg) 
{
  return atoi(arg) ;
}


////////////////////////////////////////////////////////////////////////////////
/// CINT constructor interface, return constructor string argument #idx as Double_t

Double_t RooFactoryWSTool::asDOUBLE(const char* arg) 
{
  return atof(arg) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Register foreign special objects in factory

void RooFactoryWSTool::registerSpecial(const char* typeName, RooFactoryWSTool::IFace* iface) 
{
  hooks()[typeName] = iface ;
}



////////////////////////////////////////////////////////////////////////////////

std::map<std::string,RooFactoryWSTool::IFace*>& RooFactoryWSTool::hooks() 
{
  if (_hooks) return *_hooks ;
  _hooks = new map<string,IFace*> ;
  return *_hooks ;
}



////////////////////////////////////////////////////////////////////////////////
/// Concatenate list of args into comma separated string

std::string RooFactoryWSTool::SpecialsIFace::create(RooFactoryWSTool& ft, const char* typeName, const char* instName, std::vector<std::string> args)
{
  char pargs[BUFFER_SIZE] ;
  pargs[0] = 0 ;
  vector<string>::iterator iter = args.begin() ;
  vector<string> pargv ;
  while(iter!=args.end()) {
    if (strlen(pargs)>0) strlcat(pargs,",",BUFFER_SIZE) ;
    string tmp = ft.processExpression(iter->c_str()) ;
    strlcat(pargs,tmp.c_str(),BUFFER_SIZE) ;
    pargv.push_back(tmp) ;
    ++iter ;
  }

  // Handling of special operator pdf class names
  string cl(typeName) ;
  if (cl=="SUM") {

    // SUM::name[a*A,b*B,C]
    ft.add(instName,pargs,kFALSE) ;

  } else if (cl=="RSUM") {

    // RSUM::name[a*A,b*B,C]
    ft.add(instName,pargs,kTRUE) ;

  } else if (cl=="ASUM") {

    // ASUM::name[a*A,b*B,C]
    ft.amplAdd(instName,pargs) ;

  } else if (cl=="PROD") {

    // PROD::name[A,B,C]
    ft.prod(instName,pargs) ;

  } else if (cl=="SIMUL") {

    // PROD::name[cat,state=Pdf,...]
    if (pargv.size()>1) {
      ft.simul(instName,pargv[0].c_str(),strchr(pargs,',')+1) ;
    } else {
      throw string(Form("Need at least two arguments in call to SIMUL::%s, have %d: %s",instName,(Int_t)pargv.size(),pargs)) ;
    }

  } else if (cl=="EXPR") {
    
    // EXPR::name['expr',var,var,...]
    if (args.size()<=2) {
      ft.createArg("RooGenericPdf",instName,pargs) ;
    } else {      
      char genargs[BUFFER_SIZE] ;
      strlcpy(genargs,args[0].c_str(),BUFFER_SIZE) ;      
      strlcat(genargs,",{",BUFFER_SIZE) ;
      for (UInt_t i=1 ; i<args.size() ; i++) {
	if (i!=1) strlcat(genargs,",",BUFFER_SIZE) ;
	strlcat(genargs,args[i].c_str(),BUFFER_SIZE) ;
      }
      strlcat(genargs,"}",BUFFER_SIZE) ;
      ft.createArg("RooGenericPdf",instName,genargs) ;
    }
  
  } else if (cl=="FCONV") {

    // FCONV::name[var,pdf1,pdf2]
    ft.createArg("RooFFTConvPdf",instName,pargs) ;
  
  } else if (cl=="NCONV") {

    // NCONV::name[var,pdf1,pdf2]
    ft.createArg("RooNumConvPdf",instName,pargs) ;
  
  } else if (cl=="sum") {
    
    // sum::name[a,b,c]
    ft.addfunc(instName,pargs) ;

  } else if (cl=="prod") {

    // prod::name[a,b,c]
    ft.prodfunc(instName,pargs) ;
    
  } else if (cl=="expr") {

    // expr::name['expr',var,var,...]
    if (args.size()<=2) {
      ft.createArg("RooFormulaVar",instName,pargs) ;
    } else {      
      char genargs[BUFFER_SIZE] ;
      strlcpy(genargs,args[0].c_str(),BUFFER_SIZE) ;      
      strlcat(genargs,",{",BUFFER_SIZE) ;
      for (UInt_t i=1 ; i<args.size() ; i++) {
	if (i!=1) strlcat(genargs,",",BUFFER_SIZE) ;
	strlcat(genargs,args[i].c_str(),BUFFER_SIZE) ;
      }
      strlcat(genargs,"}",BUFFER_SIZE) ;
      ft.createArg("RooFormulaVar",instName,genargs) ;
    }

  } else if (cl=="nconv") {

    // nconv::name[var,pdf1,pdf2]
    ft.createArg("RooNumConvolution",instName,pargs) ;
  
  } else if (cl=="nll") {

    // nll::name[pdf,data]
    RooNLLVar nll(instName,instName,ft.asPDF(pargv[0].c_str()),ft.asDATA(pargv[1].c_str())) ;
    if (ft.ws().import(nll,Silence())) ft.logError() ;
    
  } else if (cl=="chi2") {

    // chi2::name[pdf,data]
    RooChi2Var nll(instName,instName,ft.asPDF(pargv[0].c_str()),ft.asDHIST(pargv[1].c_str())) ;
    if (ft.ws().import(nll,Silence())) ft.logError() ;
    
  } else if (cl=="profile") {

    // profile::name[func,vars]
    ft.createArg("RooProfileLL",instName,pargs) ;

  } else if (cl=="dataobs") {

    // dataobs::name[dset,func]
    RooAbsArg* funcClone = static_cast<RooAbsArg*>(ft.asARG(pargv[1].c_str()).clone(instName)) ;
    RooAbsArg* arg = ft.asDSET(pargv[0].c_str()).addColumn(*funcClone) ;
    if (!ft.ws().fundArg(arg->GetName())) {
      if (ft.ws().import(*arg,Silence())) ft.logError() ;
    }
    delete funcClone ;

  } else if (cl=="int") {    
    
    // int::name[func,intobs]
    // int::name[func,intobs|range]
    // int::name[func,intobs,normobs]
    // int::name[func,intobs|range,normobs]

    if (pargv.size()<2 || pargv.size()>3) {
      throw string(Form("int::%s, requires 2 or 3 arguments, have %d arguments",instName,(Int_t)pargv.size())) ;
    }

    RooAbsReal& func = ft.asFUNC(pargv[0].c_str()) ;

    char buf[256] ;
    strlcpy(buf,pargv[1].c_str(),256) ;
    char* save ;
    const char* intobs = R__STRTOK_R(buf,"|",&save) ;
    if (!intobs) intobs="" ;

    const char* range = R__STRTOK_R(0,"",&save) ;
    if (!range) range="" ;

    RooAbsReal* integral = 0 ;
    if (pargv.size()==2) {
      if (range && strlen(range)) {
	integral = func.createIntegral(ft.asSET(intobs),Range(range)) ;
      } else {
	integral = func.createIntegral(ft.asSET(intobs)) ;
      }
    } else {
      if (range && strlen(range)) {
	integral = func.createIntegral(ft.asSET(intobs),Range(range),NormSet(ft.asSET(pargv[2].c_str()))) ;
      } else {
	integral = func.createIntegral(ft.asSET(intobs),NormSet(ft.asSET(pargv[2].c_str()))) ;
      }
    }

    integral->SetName(instName) ;
    if (ft.ws().import(*integral,Silence())) ft.logError() ;
    
  } else if (cl=="deriv") {
    
    // derive::name[func,obs,order]

    if (pargv.size()<2 || pargv.size()>3) {
      throw string(Form("deriv::%s, requires 2 or 3 arguments, have %d arguments",instName,(Int_t)pargv.size())) ;
    }

    RooAbsReal& func = ft.asFUNC(pargv[0].c_str()) ;

    RooAbsReal* derivative(0) ;
    if (pargv.size()==2) {
      derivative = func.derivative(ft.asVAR(pargv[1].c_str()),1) ;
    } else {
      derivative = func.derivative(ft.asVAR(pargv[1].c_str()),ft.asINT(pargv[2].c_str())) ;
    }

    derivative->SetName(instName) ;
    if (ft.ws().import(*derivative,Silence())) ft.logError() ;

  } else if (cl=="cdf") {

    // cdf::name[pdf,obs,extranormobs]

    if (pargv.size()<2 || pargv.size()>3) {
      throw string(Form("cdf::%s, requires 2 or 3 arguments, have %d arguments",instName,(Int_t)pargv.size())) ;
    }

    RooAbsPdf& pdf = ft.asPDF(pargv[0].c_str()) ;

    RooAbsReal* cdf(0) ;
    if (pargv.size()==2) {
      cdf = pdf.createCdf(ft.asSET(pargv[1].c_str())) ;
    } else {
      cdf = pdf.createCdf(ft.asSET(pargv[1].c_str()),ft.asSET(pargv[2].c_str())) ;
    }

    cdf->SetName(instName) ;
    if (ft.ws().import(*cdf,Silence())) ft.logError() ;


  } else if (cl=="PROJ") {

    // PROJ::name(pdf,intobs)
    if (pargv.size()!=2) {
      throw string(Form("PROJ::%s, requires 2 arguments, have %d arguments",instName,(Int_t)pargv.size())) ;
    }

    RooAbsPdf& pdf = ft.asPDF(pargv[0].c_str()) ;
    RooAbsPdf* projection = pdf.createProjection(ft.asSET(pargv[1].c_str())) ;
    projection->SetName(instName) ;

    if (ft.ws().import(*projection,Silence())) ft.logError() ;

  } else if (cl=="set") {

    // set::name(arg,arg,...)
    if (ft.ws().defineSet(instName,pargs)) {
      ft.logError() ;
      return string(instName) ;
    }
    
  } else {

    throw string(Form("RooFactoryWSTool::SpecialsIFace::create() ERROR: Unknown meta-type %s",typeName)) ;

  }
  return string(instName) ;    
}


RooFactoryWSTool* RooFactoryWSTool::of() 
{
  return _of ;
}

