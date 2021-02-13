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
#include "TString.h"
#include "TList.h"
#include "RooCmdArg.h"
#include "RooArgSet.h"
#include <string>

class RooCmdConfig : public TObject {
public:

  RooCmdConfig(const char* methodName);
  RooCmdConfig(const RooCmdConfig& other) ;
  ~RooCmdConfig();

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
  void defineMutex(const char* argName1, const char* argName2) ;
  void defineMutex(const char* argName1, const char* argName2, const char* argName3) ;
  void defineMutex(const char* argName1, const char* argName2, const char* argName3, const char* argName4) ;
  void defineMutex(const char* argName1, const char* argName2, const char* argName3, const char* argName4, const char* argName5) ;
  void defineRequiredArgs(const char* argName1, const char* argName2=0,
			  const char* argName3=0, const char* argName4=0,
			  const char* argName5=0, const char* argName6=0,
			  const char* argName7=0, const char* argName8=0) ;

  Bool_t defineInt(const char* name, const char* argName, Int_t intNum, Int_t defValue=0) ;
  Bool_t defineDouble(const char* name, const char* argName, Int_t doubleNum, Double_t defValue=0.) ;
  Bool_t defineString(const char* name, const char* argName, Int_t stringNum, const char* defValue="",Bool_t appendMode=kFALSE) ;
  Bool_t defineObject(const char* name, const char* argName, Int_t setNum, const TObject* obj=0, Bool_t isArray=kFALSE) ;
  Bool_t defineSet(const char* name, const char* argName, Int_t setNum, const RooArgSet* set=0) ;

  Bool_t process(const RooCmdArg& arg) ;
  Bool_t process(const RooCmdArg& arg1, const RooCmdArg& arg2, 
                 const RooCmdArg& arg3=RooCmdArg::none(), const RooCmdArg& arg4=RooCmdArg::none(), 
                 const RooCmdArg& arg5=RooCmdArg::none(), const RooCmdArg& arg6=RooCmdArg::none(), 
                 const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) ;
  Bool_t process(const RooLinkedList& argList) ;
  /// Process several RooCmdArg using iterators.
  template<typename It_t>
  bool process(It_t begin, It_t end) {
    bool result = false;
    for (auto it = begin; it != end; ++it) {
      result |= process(*it);
    }
    return result;
  }

  Int_t getInt(const char* name, Int_t defaultValue=0) ;
  Double_t getDouble(const char* name, Double_t defaultValue=0) ;
  const char* getString(const char* name, const char* defaultValue="",Bool_t convEmptyToNull=kFALSE) ;
  TObject* getObject(const char* name, TObject* obj=0) ;
  RooArgSet* getSet(const char* name, RooArgSet* set=0) ;
  const RooLinkedList& getObjectList(const char* name) ;

  Bool_t ok(Bool_t verbose) const ;

  const char* missingArgs() const ;

  RooLinkedList filterCmdList(RooLinkedList& cmdInList, const char* cmdNameList, Bool_t removeFromInList=kTRUE) ;
  void stripCmdList(RooLinkedList& cmdList, const char* cmdsToPurge) ;
  Bool_t hasProcessed(const char* cmdName) const ;

  void print() ;


  static Int_t decodeIntOnTheFly(const char* callerID, const char* cmdArgName, Int_t intIdx, Int_t defVal, const RooCmdArg& arg1, 
				 const RooCmdArg& arg2=RooCmdArg(), const RooCmdArg& arg3=RooCmdArg(), const RooCmdArg& arg4=RooCmdArg(),
				 const RooCmdArg& arg5=RooCmdArg(), const RooCmdArg& arg6=RooCmdArg(), const RooCmdArg& arg7=RooCmdArg(),
				 const RooCmdArg& arg8=RooCmdArg(), const RooCmdArg& arg9=RooCmdArg()) ;

  static std::string decodeStringOnTheFly(const char* callerID, const char* cmdArgName, Int_t intIdx, const char* defVal, const RooCmdArg& arg1,
					 const RooCmdArg& arg2=RooCmdArg(), const RooCmdArg& arg3=RooCmdArg(), const RooCmdArg& arg4=RooCmdArg(),
					 const RooCmdArg& arg5=RooCmdArg(), const RooCmdArg& arg6=RooCmdArg(), const RooCmdArg& arg7=RooCmdArg(),
					 const RooCmdArg& arg8=RooCmdArg(), const RooCmdArg& arg9=RooCmdArg()) ;

  static TObject* decodeObjOnTheFly(const char* callerID, const char* cmdArgName, Int_t objIdx, TObject* defVal, const RooCmdArg& arg1, 
				     const RooCmdArg& arg2=RooCmdArg(), const RooCmdArg& arg3=RooCmdArg(), const RooCmdArg& arg4=RooCmdArg(),
				     const RooCmdArg& arg5=RooCmdArg(), const RooCmdArg& arg6=RooCmdArg(), const RooCmdArg& arg7=RooCmdArg(),
				     const RooCmdArg& arg8=RooCmdArg(), const RooCmdArg& arg9=RooCmdArg()) ;

  static double decodeDoubleOnTheFly(const char* callerID, const char* cmdArgName, int idx, double defVal,
      std::initializer_list<std::reference_wrapper<const RooCmdArg>> args);

protected:

  TString _name ;
  
  Bool_t _verbose ;
  Bool_t _error ;
  Bool_t _allowUndefined ;

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

  TIterator* _iIter ; // Iterator over integer list
  TIterator* _dIter ; // Iterator over double list
  TIterator* _sIter ; // Iterator over string list
  TIterator* _oIter ; // Iterator over object list
  TIterator* _cIter ; // Iterator over RooArgSet list
  TIterator* _rIter ; // Iterator over required cmd list
  TIterator* _fIter ; // Iterator over forbidden cmd list
  TIterator* _mIter ; // Iterator over mutex list
  TIterator* _yIter ; // Iterator over dependency list
  TIterator* _pIter ; // Iterator over processed cmd list

  ClassDef(RooCmdConfig,0) // Configurable parse of RooCmdArg objects
};

#endif


