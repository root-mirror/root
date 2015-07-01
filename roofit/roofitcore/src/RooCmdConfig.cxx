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
// 
// BEGIN_HTML
// Class RooCmdConfig is a configurable parser for RooCmdArg named
// arguments. It maps the contents of named arguments named to integers,
// doubles, strings and TObjects that can be retrieved after processing
// a set of RooCmdArgs. The parser also has options to enforce syntax
// rules such as (conditionally) required arguments, mutually exclusive
// arguments and dependencies between arguments
// END_HTML
//

#include "RooFit.h"

#include "RooCmdConfig.h"
#include "RooInt.h"
#include "RooDouble.h"
#include "RooArgSet.h"
#include "RooStringVar.h"
#include "RooTObjWrap.h"
#include "RooAbsData.h"
#include "TObjString.h"
#include "RooMsgService.h"

#include "Riostream.h"


using namespace std;

ClassImp(RooCmdConfig) 
  ;



////////////////////////////////////////////////////////////////////////////////
/// Constructor taking descriptive name of owner/user which
/// is used as prefix for any warning or error messages
/// generated by this parser

RooCmdConfig::RooCmdConfig(const char* methodName) :
  TObject(),
  _name(methodName)
{
  _verbose = kFALSE ;
  _error = kFALSE ;
  _allowUndefined = kFALSE ;

  _iIter = _iList.MakeIterator() ;
  _dIter = _dList.MakeIterator() ;
  _sIter = _sList.MakeIterator() ;
  _oIter = _oList.MakeIterator() ;
  _cIter = _cList.MakeIterator() ;

  _rIter = _rList.MakeIterator() ;
  _fIter = _fList.MakeIterator() ;
  _mIter = _mList.MakeIterator() ;
  _yIter = _yList.MakeIterator() ;
  _pIter = _pList.MakeIterator() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooCmdConfig::RooCmdConfig(const RooCmdConfig& other)  : TObject(other)
{
  _name   = other._name ;
  _verbose = other._verbose ;
  _error = other._error ;
  _allowUndefined = other._allowUndefined ;

  _iIter = _iList.MakeIterator() ;
  _dIter = _dList.MakeIterator() ;
  _sIter = _sList.MakeIterator() ;
  _oIter = _oList.MakeIterator() ;
  _cIter = _cList.MakeIterator() ;
  _rIter = _rList.MakeIterator() ;
  _fIter = _fList.MakeIterator() ;
  _mIter = _mList.MakeIterator() ;
  _yIter = _yList.MakeIterator() ;
  _pIter = _pList.MakeIterator() ;

  other._iIter->Reset() ;
  RooInt* ri ;
  while((ri=(RooInt*)other._iIter->Next())) {
    _iList.Add(ri->Clone()) ;
  }

  other._dIter->Reset() ;
  RooDouble* rd ;
  while((rd=(RooDouble*)other._dIter->Next())) {
    _dList.Add(rd->Clone()) ;
  }

  other._sIter->Reset() ;
  RooStringVar* rs ;
  while((rs=(RooStringVar*)other._sIter->Next())) {
    _sList.Add(rs->Clone()) ;
  }

  other._oIter->Reset() ;
  RooTObjWrap* os ;
  while((os=(RooTObjWrap*)other._oIter->Next())) {
    _oList.Add(os->Clone()) ;
  }

  other._cIter->Reset() ;
  RooTObjWrap* cs ;
  while((cs=(RooTObjWrap*)other._cIter->Next())) {
    _cList.Add(cs->Clone()) ;
  }

  other._rIter->Reset() ;
  TObjString* rr ;
  while((rr=(TObjString*)other._rIter->Next())) {
    _rList.Add(rr->Clone()) ;
  }

  other._fIter->Reset() ;
  TObjString* ff ;
  while((ff=(TObjString*)other._fIter->Next())) {
    _fList.Add(ff->Clone()) ;
  }

  other._mIter->Reset() ;
  TObjString* mm ;
  while((mm=(TObjString*)other._mIter->Next())) {
    _mList.Add(mm->Clone()) ;
  }

  other._yIter->Reset() ;
  TObjString* yy ;
  while((yy=(TObjString*)other._yIter->Next())) {
    _yList.Add(yy->Clone()) ;
  }

  other._pIter->Reset() ;
  TObjString* pp ;
  while((pp=(TObjString*)other._pIter->Next())) {
    _pList.Add(pp->Clone()) ;
  }

}



////////////////////////////////////////////////////////////////////////////////
/// Destructor 

RooCmdConfig::~RooCmdConfig()
{
  delete _iIter ;
  delete _dIter ;
  delete _sIter ;
  delete _oIter ;
  delete _cIter ;
  delete _rIter ;
  delete _fIter ;
  delete _mIter ;
  delete _yIter ;
  delete _pIter ;

  _iList.Delete() ;
  _dList.Delete() ;
  _sList.Delete() ;
  _cList.Delete() ;
  _oList.Delete() ;
  _rList.Delete() ;
  _fList.Delete() ;
  _mList.Delete() ;
  _yList.Delete() ;
  _pList.Delete() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Add condition that any of listed arguments must be processed
/// for parsing to be declared successful

void RooCmdConfig::defineRequiredArgs(const char* argName1, const char* argName2,
				      const char* argName3, const char* argName4,
				      const char* argName5, const char* argName6,
				      const char* argName7, const char* argName8) 
{
  if (argName1) _rList.Add(new TObjString(argName1)) ;
  if (argName2) _rList.Add(new TObjString(argName2)) ;
  if (argName3) _rList.Add(new TObjString(argName3)) ;
  if (argName4) _rList.Add(new TObjString(argName4)) ;
  if (argName5) _rList.Add(new TObjString(argName5)) ;
  if (argName6) _rList.Add(new TObjString(argName6)) ;
  if (argName7) _rList.Add(new TObjString(argName7)) ;
  if (argName8) _rList.Add(new TObjString(argName8)) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return string with names of arguments that were required, but not
/// processed

const char* RooCmdConfig::missingArgs() const 
{
  static TString ret ;
  ret="" ;

  _rIter->Reset() ;
  TObjString* s ;
  Bool_t first(kTRUE) ;
  while((s=(TObjString*)_rIter->Next())) {
    if (first) {
      first=kFALSE ;
    } else {
      ret.Append(", ") ;
    }
    ret.Append(s->String()) ;
  }

  return ret.Length() ? ret.Data() : 0 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Define that processing argument name refArgName requires processing
/// of argument named neededArgName to successfully complete parsing

void RooCmdConfig::defineDependency(const char* refArgName, const char* neededArgName) 
{
  TNamed* dep = new TNamed(refArgName,neededArgName) ;
  _yList.Add(dep) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Define arguments named argName1 and argName2 mutually exclusive

void RooCmdConfig::defineMutex(const char* argName1, const char* argName2) 
{
  TNamed* mutex1 = new TNamed(argName1,argName2) ;
  TNamed* mutex2 = new TNamed(argName2,argName1) ;
  _mList.Add(mutex1) ;
  _mList.Add(mutex2) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Define arguments named argName1,argName2 and argName3 mutually exclusive

void RooCmdConfig::defineMutex(const char* argName1, const char* argName2, const char* argName3) 
{
  defineMutex(argName1,argName2) ;
  defineMutex(argName1,argName3) ;
  defineMutex(argName2,argName3) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Define arguments named argName1,argName2,argName3 and argName4 mutually exclusive

void RooCmdConfig::defineMutex(const char* argName1, const char* argName2, const char* argName3, const char* argName4) 
{
  defineMutex(argName1,argName2) ;
  defineMutex(argName1,argName3) ;
  defineMutex(argName1,argName4) ;
  defineMutex(argName2,argName3) ;
  defineMutex(argName2,argName4) ;
  defineMutex(argName3,argName4) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Define arguments named argName1,argName2,argName3 and argName4 mutually exclusive

void RooCmdConfig::defineMutex(const char* argName1, const char* argName2, const char* argName3, const char* argName4, const char* argName5) 
{
  defineMutex(argName1,argName2) ;
  defineMutex(argName1,argName3) ;
  defineMutex(argName1,argName4) ;
  defineMutex(argName1,argName4) ;
  defineMutex(argName2,argName3) ;
  defineMutex(argName2,argName4) ;
  defineMutex(argName2,argName4) ;
  defineMutex(argName3,argName4) ;
  defineMutex(argName3,argName5) ;
  defineMutex(argName4,argName5) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Define integer property name 'name' mapped to integer in slot 'intNum' in RooCmdArg with name argName
/// Define default value for this int property to be defVal in case named argument is not processed

Bool_t RooCmdConfig::defineInt(const char* name, const char* argName, Int_t intNum, Int_t defVal)
{
  if (_iList.FindObject(name)) {
    coutE(InputArguments) << "RooCmdConfig::defintInt: name '" << name << "' already defined" << endl ;
    return kTRUE ;
  }

  RooInt* ri = new RooInt(defVal) ;
  ri->SetName(name) ;
  ri->SetTitle(argName) ;
  ri->SetUniqueID(intNum) ;
  
  _iList.Add(ri) ;
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Define Double_t property name 'name' mapped to Double_t in slot 'doubleNum' in RooCmdArg with name argName
/// Define default value for this Double_t property to be defVal in case named argument is not processed

Bool_t RooCmdConfig::defineDouble(const char* name, const char* argName, Int_t doubleNum, Double_t defVal) 
{
  if (_dList.FindObject(name)) {
    coutE(InputArguments) << "RooCmdConfig::defineDouble: name '" << name << "' already defined" << endl ;
    return kTRUE ;
  }

  RooDouble* rd = new RooDouble(defVal) ;
  rd->SetName(name) ;
  rd->SetTitle(argName) ;
  rd->SetUniqueID(doubleNum) ;
  
  _dList.Add(rd) ;
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Define Double_t property name 'name' mapped to Double_t in slot 'stringNum' in RooCmdArg with name argName
/// Define default value for this Double_t property to be defVal in case named argument is not processed
/// If appendMode is true, values found in multiple matching RooCmdArg arguments will be concatenated
/// in the output string. If it is false, only the value of the last processed instance is retained

Bool_t RooCmdConfig::defineString(const char* name, const char* argName, Int_t stringNum, const char* defVal, Bool_t appendMode) 
{
  if (_sList.FindObject(name)) {
    coutE(InputArguments) << "RooCmdConfig::defineString: name '" << name << "' already defined" << endl ;
    return kTRUE ;
  }

  RooStringVar* rs = new RooStringVar(name,argName,defVal,64000) ;
  if (appendMode) {
    rs->setAttribute("RooCmdConfig::AppendMode") ;
  }
  rs->SetUniqueID(stringNum) ;
  
  _sList.Add(rs) ;
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Define TObject property name 'name' mapped to object in slot 'setNum' in RooCmdArg with name argName
/// Define default value for this TObject property to be defVal in case named argument is not processed.
/// If isArray is true, an array of TObjects is harvested in case multiple matching named arguments are processed.
/// If isArray is false, only the TObject in the last processed named argument is retained

Bool_t RooCmdConfig::defineObject(const char* name, const char* argName, Int_t setNum, const TObject* defVal, Bool_t isArray) 
{

  if (_oList.FindObject(name)) {
    coutE(InputArguments) << "RooCmdConfig::defineObject: name '" << name << "' already defined" << endl ;
    return kTRUE ;
  }

  RooTObjWrap* os = new RooTObjWrap((TObject*)defVal,isArray) ;
  os->SetName(name) ;
  os->SetTitle(argName) ;
  os->SetUniqueID(setNum) ;
  
  _oList.Add(os) ;
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Define TObject property name 'name' mapped to object in slot 'setNum' in RooCmdArg with name argName
/// Define default value for this TObject property to be defVal in case named argument is not processed.
/// If isArray is true, an array of TObjects is harvested in case multiple matching named arguments are processed.
/// If isArray is false, only the TObject in the last processed named argument is retained

Bool_t RooCmdConfig::defineSet(const char* name, const char* argName, Int_t setNum, const RooArgSet* defVal) 
{

  if (_cList.FindObject(name)) {
    coutE(InputArguments) << "RooCmdConfig::defineObject: name '" << name << "' already defined" << endl ;
    return kTRUE ;
  }

  RooTObjWrap* cs = new RooTObjWrap((TObject*)defVal) ;
  cs->SetName(name) ;
  cs->SetTitle(argName) ;
  cs->SetUniqueID(setNum) ;
  
  _cList.Add(cs) ;
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print configuration of parser

void RooCmdConfig::print()
{
  // Find registered integer fields for this opcode 
  _iIter->Reset() ;
  RooInt* ri ;
  while((ri=(RooInt*)_iIter->Next())) {
    cout << ri->GetName() << "[Int_t] = " << *ri << endl ;
  }

  // Find registered double fields for this opcode 
  _dIter->Reset() ;
  RooDouble* rd ;
  while((rd=(RooDouble*)_dIter->Next())) {
    cout << rd->GetName() << "[Double_t] = " << *rd << endl ;
  }

  // Find registered string fields for this opcode 
  _sIter->Reset() ;
  RooStringVar* rs ;
  while((rs=(RooStringVar*)_sIter->Next())) {
    cout << rs->GetName() << "[string] = \"" << rs->getVal() << "\"" << endl ;
  }

  // Find registered argset fields for this opcode 
  _oIter->Reset() ;
  RooTObjWrap* ro ;
  while((ro=(RooTObjWrap*)_oIter->Next())) {
    cout << ro->GetName() << "[TObject] = " ; 
    if (ro->obj()) {
      cout << ro->obj()->GetName() << endl ;
    } else {

      cout << "(null)" << endl ;
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Process given list with RooCmdArgs

Bool_t RooCmdConfig::process(const RooLinkedList& argList) 
{
  Bool_t ret(kFALSE) ;
  TIterator* iter = argList.MakeIterator() ;
  RooCmdArg* arg ;
  while((arg=(RooCmdArg*)iter->Next())) {
    ret |= process(*arg) ;
  }
  delete iter ;
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Process given RooCmdArgs

Bool_t RooCmdConfig::process(const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3, const RooCmdArg& arg4,
			     const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8) 
{
  Bool_t ret(kFALSE) ;
  ret |= process(arg1) ;
  ret |= process(arg2) ;
  ret |= process(arg3) ;
  ret |= process(arg4) ;
  ret |= process(arg5) ;
  ret |= process(arg6) ;
  ret |= process(arg7) ;
  ret |= process(arg8) ;
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Process given RooCmdArg

Bool_t RooCmdConfig::process(const RooCmdArg& arg) 
{
  // Retrive command code
  const char* opc = arg.opcode() ;

  // Ignore empty commands
  if (!opc) return kFALSE ;

  // Check if not forbidden
  if (_fList.FindObject(opc)) {
    coutE(InputArguments) << _name << " ERROR: argument " << opc << " not allowed in this context" << endl ;
    _error = kTRUE ;
    return kTRUE ;
  }

  // Check if this code generates any dependencies
  TObject* dep = _yList.FindObject(opc) ;
  if (dep) {
    // Dependent command found, add to required list if not already processed
    if (!_pList.FindObject(dep->GetTitle())) {
      _rList.Add(new TObjString(dep->GetTitle())) ;
      if (_verbose) {
	cout << "RooCmdConfig::process: " << opc << " has unprocessed dependent " << dep->GetTitle() 
	     << ", adding to required list" << endl ;
      }
    } else {
      if (_verbose) {
	cout << "RooCmdConfig::process: " << opc << " dependent " << dep->GetTitle() << " is already processed" << endl ;
      }
    }
  }

  // Check for mutexes
  TObject * mutex = _mList.FindObject(opc) ;
  if (mutex) {
    if (_verbose) {
      cout << "RooCmdConfig::process: " << opc << " excludes " << mutex->GetTitle() 
	   << ", adding to forbidden list" << endl ;
    }    
    _fList.Add(new TObjString(mutex->GetTitle())) ;
  }


  Bool_t anyField(kFALSE) ;

  // Find registered integer fields for this opcode 
  _iIter->Reset() ;
  RooInt* ri ;
  while((ri=(RooInt*)_iIter->Next())) {
    if (!TString(opc).CompareTo(ri->GetTitle())) {
      *ri = arg.getInt(ri->GetUniqueID()) ;
      anyField = kTRUE ;
      if (_verbose) {
	cout << "RooCmdConfig::process " << ri->GetName() << "[Int_t]" << " set to " << *ri << endl ;
      }
    }
  }

  // Find registered double fields for this opcode 
  _dIter->Reset() ;
  RooDouble* rd ;
  while((rd=(RooDouble*)_dIter->Next())) {
    if (!TString(opc).CompareTo(rd->GetTitle())) {
      *rd = arg.getDouble(rd->GetUniqueID()) ;
      anyField = kTRUE ;
      if (_verbose) {
	cout << "RooCmdConfig::process " << rd->GetName() << "[Double_t]" << " set to " << *rd << endl ;
      }
    }
  }

  // Find registered string fields for this opcode 
  _sIter->Reset() ;
  RooStringVar* rs ;
  while((rs=(RooStringVar*)_sIter->Next())) {
    if (!TString(opc).CompareTo(rs->GetTitle())) {
      
      const char* oldStr = rs->getVal() ;

      if (oldStr && strlen(oldStr)>0 && rs->getAttribute("RooCmdConfig::AppendMode")) {
	rs->setVal(Form("%s,%s",rs->getVal(),arg.getString(rs->GetUniqueID()))) ;
      } else {
	rs->setVal(arg.getString(rs->GetUniqueID())) ;
      }
      anyField = kTRUE ;
      if (_verbose) {
	cout << "RooCmdConfig::process " << rs->GetName() << "[string]" << " set to " << rs->getVal() << endl ;
      }
    }
  }

  // Find registered TObject fields for this opcode 
  _oIter->Reset() ;
  RooTObjWrap* os ;
  while((os=(RooTObjWrap*)_oIter->Next())) {
    if (!TString(opc).CompareTo(os->GetTitle())) {
      os->setObj((TObject*)arg.getObject(os->GetUniqueID())) ;
      anyField = kTRUE ;
      if (_verbose) {
	cout << "RooCmdConfig::process " << os->GetName() << "[TObject]" << " set to " ;
	if (os->obj()) {
	  cout << os->obj()->GetName() << endl ;
	} else {
	  cout << "(null)" << endl ;
	}
      }
    }
  }

  // Find registered RooArgSet fields for this opcode 
  _cIter->Reset() ;
  RooTObjWrap* cs ;
  while((cs=(RooTObjWrap*)_cIter->Next())) {
    if (!TString(opc).CompareTo(cs->GetTitle())) {
      cs->setObj((TObject*)arg.getSet(cs->GetUniqueID())) ;
      anyField = kTRUE ;
      if (_verbose) {
	cout << "RooCmdConfig::process " << cs->GetName() << "[RooArgSet]" << " set to " ;
	if (cs->obj()) {
	  cout << cs->obj()->GetName() << endl ;
	} else {
	  cout << "(null)" << endl ;
	}
      }
    }
  }

  Bool_t multiArg = !TString("MultiArg").CompareTo(opc) ;

  if (!anyField && !_allowUndefined && !multiArg) {
    coutE(InputArguments) << _name << " ERROR: unrecognized command: " << opc << endl ;
  }


  // Remove command from required-args list (if it was there)
  TObject* obj = _rList.FindObject(opc) ;
  if (obj) {
    _rList.Remove(obj) ;
  }

  // Add command the processed list
  TNamed *pcmd = new TNamed(opc,opc) ;
  _pList.Add(pcmd) ;

  Bool_t depRet = kFALSE ;
  if (arg._procSubArgs) {
    for (Int_t ia=0 ; ia<arg._argList.GetSize() ; ia++) {
      RooCmdArg* subArg = static_cast<RooCmdArg*>(arg._argList.At(ia)) ;
      if (strlen(subArg->GetName())>0) {
	RooCmdArg subArgCopy(*subArg) ;
	if (arg._prefixSubArgs) {
	  subArgCopy.SetName(Form("%s::%s",arg.GetName(),subArg->GetName())) ;
	}
	depRet |= process(subArgCopy) ;
      }
    }
  }

  return ((anyField||_allowUndefined)?kFALSE:kTRUE)||depRet ;
}
  


////////////////////////////////////////////////////////////////////////////////
/// Return true if RooCmdArg with name 'cmdName' has been processed

Bool_t RooCmdConfig::hasProcessed(const char* cmdName) const 
{
  return _pList.FindObject(cmdName) ? kTRUE : kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return integer property registered with name 'name'. If no
/// property is registered, return defVal

Int_t RooCmdConfig::getInt(const char* name, Int_t defVal) 
{
  RooInt* ri = (RooInt*) _iList.FindObject(name) ;
  return ri ? (Int_t)(*ri) : defVal ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return Double_t property registered with name 'name'. If no
/// property is registered, return defVal

Double_t RooCmdConfig::getDouble(const char* name, Double_t defVal) 
{
  RooDouble* rd = (RooDouble*) _dList.FindObject(name) ;
  return rd ? (Double_t)(*rd) : defVal ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return string property registered with name 'name'. If no
/// property is registered, return defVal. If convEmptyToNull
/// is true, empty string will be returned as null pointers

const char* RooCmdConfig::getString(const char* name, const char* defVal, Bool_t convEmptyToNull) 
{
  RooStringVar* rs = (RooStringVar*) _sList.FindObject(name) ;
  return rs ? ((convEmptyToNull && strlen(rs->getVal())==0) ? 0 : ((const char*)rs->getVal()) ) : defVal ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return TObject property registered with name 'name'. If no
/// property is registered, return defVal

TObject* RooCmdConfig::getObject(const char* name, TObject* defVal) 
{
  RooTObjWrap* ro = (RooTObjWrap*) _oList.FindObject(name) ;
  return ro ? ro->obj() : defVal ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return RooArgSet property registered with name 'name'. If no
/// property is registered, return defVal

RooArgSet* RooCmdConfig::getSet(const char* name, RooArgSet* defVal) 
{
  RooTObjWrap* ro = (RooTObjWrap*) _cList.FindObject(name) ;
  return ro ? ((RooArgSet*)ro->obj()) : defVal ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return list of objects registered with name 'name'

const RooLinkedList& RooCmdConfig::getObjectList(const char* name) 
{
  static RooLinkedList defaultDummy ;
  RooTObjWrap* ro = (RooTObjWrap*) _oList.FindObject(name) ;
  return ro ? ro->objList() : defaultDummy ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return true of parsing was successful

Bool_t RooCmdConfig::ok(Bool_t verbose) const 
{ 
  if (_rList.GetSize()==0 && !_error) return kTRUE ;

  if (verbose) {
    const char* margs = missingArgs() ;
    if (margs) {
      coutE(InputArguments) << _name << " ERROR: missing arguments: " << margs << endl ;
    } else {
      coutE(InputArguments) << _name << " ERROR: illegal combination of arguments and/or missing arguments" << endl ;
    }
  }
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Utility function that strips command names listed (comma separated) in cmdsToPurge from cmdList

void RooCmdConfig::stripCmdList(RooLinkedList& cmdList, const char* cmdsToPurge) 
{
  // Sanity check
  if (!cmdsToPurge) return ;
  
  // Copy command list for parsing
  char buf[1024] ;
  strlcpy(buf,cmdsToPurge,1024) ;
  
  char* name = strtok(buf,",") ;
  while(name) {
    TObject* cmd = cmdList.FindObject(name) ;
    if (cmd) cmdList.Remove(cmd) ;
    name = strtok(0,",") ;
  }

}



////////////////////////////////////////////////////////////////////////////////
/// Utility function to filter commands listed in cmdNameList from cmdInList. Filtered arguments are put in the returned list.
/// If removeFromInList is true then these commands are removed from the input list

RooLinkedList RooCmdConfig::filterCmdList(RooLinkedList& cmdInList, const char* cmdNameList, Bool_t removeFromInList) 
{
  RooLinkedList filterList ;
  if (!cmdNameList) return filterList ;

  // Copy command list for parsing
  char buf[1024] ;
  strlcpy(buf,cmdNameList,1024) ;
  
  char* name = strtok(buf,",") ;
  while(name) {
    TObject* cmd = cmdInList.FindObject(name) ;
    if (cmd) {
      if (removeFromInList) {
	cmdInList.Remove(cmd) ;
      }
      filterList.Add(cmd) ;
    }
    name = strtok(0,",") ;
  }
  return filterList ;  
}



////////////////////////////////////////////////////////////////////////////////
/// Static decoder function allows to retrieve integer property from set of RooCmdArgs 
/// For use in base member initializers in constructors

Int_t RooCmdConfig::decodeIntOnTheFly(const char* callerID, const char* cmdArgName, Int_t intIdx, Int_t defVal, const RooCmdArg& arg1, 
				      const RooCmdArg& arg2, const RooCmdArg& arg3, const RooCmdArg& arg4,
				      const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7,
				      const RooCmdArg& arg8, const RooCmdArg& arg9) 
{
  RooCmdConfig pc(callerID) ;
  pc.allowUndefined() ;
  pc.defineInt("theInt",cmdArgName,intIdx,defVal) ;
  pc.process(arg1) ;  pc.process(arg2) ;  pc.process(arg3) ;
  pc.process(arg4) ;  pc.process(arg5) ;  pc.process(arg6) ;
  pc.process(arg7) ;  pc.process(arg8) ;  pc.process(arg9) ;
  return pc.getInt("theInt") ;
}



////////////////////////////////////////////////////////////////////////////////
/// Static decoder function allows to retrieve string property from set of RooCmdArgs 
/// For use in base member initializers in constructors

const char* RooCmdConfig::decodeStringOnTheFly(const char* callerID, const char* cmdArgName, Int_t strIdx, const char* defVal, const RooCmdArg& arg1, 
					 const RooCmdArg& arg2, const RooCmdArg& arg3, const RooCmdArg& arg4,
					 const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7,
					 const RooCmdArg& arg8, const RooCmdArg& arg9) 
{  
  static string retBuf = "" ;

  RooCmdConfig pc(callerID) ;
  pc.allowUndefined() ;
  pc.defineString("theString",cmdArgName,strIdx,defVal) ;
  pc.process(arg1) ;  pc.process(arg2) ;  pc.process(arg3) ;
  pc.process(arg4) ;  pc.process(arg5) ;  pc.process(arg6) ;
  pc.process(arg7) ;  pc.process(arg8) ;  pc.process(arg9) ;
  const char* ret =  pc.getString("theString",0,kTRUE) ;

  if (ret) {
    retBuf = ret ;
  } else {
    retBuf.clear() ;
  }
  return retBuf.c_str() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Static decoder function allows to retrieve object property from set of RooCmdArgs 
/// For use in base member initializers in constructors

TObject* RooCmdConfig::decodeObjOnTheFly(const char* callerID, const char* cmdArgName, Int_t objIdx, TObject* defVal, const RooCmdArg& arg1, 
					 const RooCmdArg& arg2, const RooCmdArg& arg3, const RooCmdArg& arg4,
					 const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7,
					 const RooCmdArg& arg8, const RooCmdArg& arg9) 
{
  RooCmdConfig pc(callerID) ;
  pc.allowUndefined() ;
  pc.defineObject("theObj",cmdArgName,objIdx,defVal) ;
  pc.process(arg1) ;  pc.process(arg2) ;  pc.process(arg3) ;
  pc.process(arg4) ;  pc.process(arg5) ;  pc.process(arg6) ;
  pc.process(arg7) ;  pc.process(arg8) ;  pc.process(arg9) ;
  return (TObject*) pc.getObject("theObj") ;
}
