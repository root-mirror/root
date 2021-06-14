/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooCmdConfig.h,v 1.12 2007/05/11 09:11:30 verkerke Exp $
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

#ifndef ROO_CMD_CONFIG
#define ROO_CMD_CONFIG

#include "TObject.h"
#include "TObjString.h"
#include "TString.h"
#include "TList.h"
#include "RooCmdArg.h"
#include "RooArgSet.h"
#include <string>

class RooCmdConfig : public TObject {
public:

  RooCmdConfig(const char* methodName);
  RooCmdConfig(const RooCmdConfig& other) ;

  void setVerbose(Bool_t flag) { 
    // If flag is true verbose messaging is activated
    _verbose = flag ; 
  }

  void allowUndefined(Bool_t flag=kTRUE) { 
    // If flag is true the processing of unrecognized RooCmdArgs
    // is not considered an error
    _allowUndefined = flag ; 
  }
  void defineDependency(const char* refArgName, const char* neededArgName) ;

  template<class... Args_t>
  void defineRequiredArgs(const char* first, Args_t && ... args);

  template<class... Args_t>
  void defineMutex(const char* head, Args_t && ... tail);
  void defineMutex(const char*) {} // to end the recursion of defineMutex()

  Bool_t defineInt(const char* name, const char* argName, Int_t intNum, Int_t defValue=0) ;
  Bool_t defineDouble(const char* name, const char* argName, Int_t doubleNum, Double_t defValue=0.) ;
  Bool_t defineString(const char* name, const char* argName, Int_t stringNum, const char* defValue="",Bool_t appendMode=kFALSE) ;
  Bool_t defineObject(const char* name, const char* argName, Int_t setNum, const TObject* obj=0, Bool_t isArray=kFALSE) ;
  Bool_t defineSet(const char* name, const char* argName, Int_t setNum, const RooArgSet* set=0) ;

  Bool_t process(const RooCmdArg& arg) ;
  template<class... Args_t>
  bool process(const RooCmdArg& arg, Args_t && ...args);
  Bool_t process(const RooLinkedList& argList) ;
  template<typename It_t>
  bool process(It_t begin, It_t end);

  Int_t getInt(const char* name, Int_t defaultValue=0) ;
  Double_t getDouble(const char* name, Double_t defaultValue=0) ;
  const char* getString(const char* name, const char* defaultValue="",Bool_t convEmptyToNull=kFALSE) ;
  TObject* getObject(const char* name, TObject* obj=0) ;
  RooArgSet* getSet(const char* name, RooArgSet* set=0) ;
  const RooLinkedList& getObjectList(const char* name) ;

  Bool_t ok(Bool_t verbose) const ;

  std::string missingArgs() const ;

  RooLinkedList filterCmdList(RooLinkedList& cmdInList, const char* cmdNameList, bool removeFromInList=true) const;
  void stripCmdList(RooLinkedList& cmdList, const char* cmdsToPurge) const;
  Bool_t hasProcessed(const char* cmdName) const ;

  void print() ;


  template<class ...Args_t>
  static Int_t decodeIntOnTheFly(
          const char* callerID, const char* cmdArgName, Int_t intIdx, Int_t defVal, Args_t && ...args);

  template<class ...Args_t>
  static std::string decodeStringOnTheFly(
          const char* callerID, const char* cmdArgName, Int_t intIdx, const char* defVal, Args_t && ...args);

  template<class ...Args_t>
  static TObject* decodeObjOnTheFly(
          const char* callerID, const char* cmdArgName, Int_t objIdx, TObject* defVal, Args_t && ...args);

  static double decodeDoubleOnTheFly(const char* callerID, const char* cmdArgName, int idx, double defVal,
      std::initializer_list<std::reference_wrapper<const RooCmdArg>> args);

protected:

  TString _name ;
  
  Bool_t _verbose = false;
  Bool_t _error = false;
  Bool_t _allowUndefined = false;

  TList _iList ; // Integer list
  TList _dList ; // Double list
  TList _sList ; // String list
  TList _oList ; // Object list
  TList _cList ; // RooArgSet list

  TList _rList ; // Required cmd list
  TList _fList ; // Forbidden cmd list
  TList _mList ; // Mutex cmd list 
  TList _yList ; // Dependency cmd list
  TList _pList ; // Processed cmd list 

  ClassDef(RooCmdConfig,1) // Configurable parse of RooCmdArg objects
};


////////////////////////////////////////////////////////////////////////////////
/// Add condition that any of listed arguments must be processed
/// for parsing to be declared successful
template<class... Args_t>
void RooCmdConfig::defineRequiredArgs(const char* first, Args_t && ... args) {
  for(auto const& arg : {first, args...}) {
      if (arg) _rList.Add(new TObjString(arg));
  }
}


//////////////////////////////////////////////////////////////////////////////////
/// Define arguments where any pair is mutually exclusive
template<class... Args_t>
void RooCmdConfig::defineMutex(const char* head, Args_t && ... tail) {
  for(auto const& item : {tail...}) {
    _mList.Add(new TNamed(head,item));
    _mList.Add(new TNamed(item,head));
  }
  defineMutex(tail...);
}


////////////////////////////////////////////////////////////////////////////////
/// Process given RooCmdArgs
template<class... Args_t>
bool RooCmdConfig::process(const RooCmdArg& arg, Args_t && ...args) {
  bool result = false;
  for(auto r : {process(arg), process(std::forward<Args_t>(args))...}) result |= r;
  return result;
}


////////////////////////////////////////////////////////////////////////////////
/// Process several RooCmdArg using iterators.
template<typename It_t>
bool RooCmdConfig::process(It_t begin, It_t end) {
  bool result = false;
  for (auto it = begin; it != end; ++it) {
    result |= process(*it);
  }
  return result;
}


////////////////////////////////////////////////////////////////////////////////
/// Static decoder function allows to retrieve integer property from set of RooCmdArgs 
/// For use in base member initializers in constructors

template<class ...Args_t>
Int_t RooCmdConfig::decodeIntOnTheFly(
        const char* callerID, const char* cmdArgName, Int_t intIdx, Int_t defVal, Args_t && ...args)
{
  RooCmdConfig pc(callerID) ;
  pc.allowUndefined() ;
  pc.defineInt("theInt",cmdArgName,intIdx,defVal) ;
  pc.process(std::forward<Args_t>(args)...);
  return pc.getInt("theInt") ;
}


////////////////////////////////////////////////////////////////////////////////
/// Static decoder function allows to retrieve string property from set of RooCmdArgs 
/// For use in base member initializers in constructors

template<class ...Args_t>
std::string RooCmdConfig::decodeStringOnTheFly(
        const char* callerID, const char* cmdArgName, Int_t strIdx, const char* defVal, Args_t && ...args)
{  
  RooCmdConfig pc(callerID) ;
  pc.allowUndefined() ;
  pc.defineString("theString",cmdArgName,strIdx,defVal) ;
  pc.process(std::forward<Args_t>(args)...);
  const char* ret =  pc.getString("theString",0,kTRUE) ;

  return ret ? ret : "";
}


////////////////////////////////////////////////////////////////////////////////
/// Static decoder function allows to retrieve object property from set of RooCmdArgs 
/// For use in base member initializers in constructors

template<class ...Args_t>
TObject* RooCmdConfig::decodeObjOnTheFly(
        const char* callerID, const char* cmdArgName, Int_t objIdx, TObject* defVal, Args_t && ...args)
{
  RooCmdConfig pc(callerID) ;
  pc.allowUndefined() ;
  pc.defineObject("theObj",cmdArgName,objIdx,defVal) ;
  pc.process(std::forward<Args_t>(args)...);
  return pc.getObject("theObj") ;
}


#endif
