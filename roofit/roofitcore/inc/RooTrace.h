/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooTrace.h,v 1.16 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_TRACE
#define ROO_TRACE

#include "RooLinkedList.h"
#include <map>
#include <string>

#define TRACE_CREATE
#define TRACE_DESTROY

class RooTrace {
public:
  RooTrace() ;
  virtual ~RooTrace() {} ;

  static void create(const TObject* obj) ;
  static void destroy(const TObject* obj) ;


  static void createSpecial(const char* name, int size) ;
  static void destroySpecial(const char* name) ;

  
  static void active(Bool_t flag) ;
  static void verbose(Bool_t flag) ;

  static void dump() ;
  static void dump(std::ostream& os, Bool_t sinceMarked=kFALSE) ;
  static void mark() ;
  
  static void callgrind_zero() ;
  static void callgrind_dump() ;

  
  static RooTrace& instance() ;

  static void printObjectCounts() ;

  static std::map<std::string,double> objectTiming;

  // TODO: for windows version, change to private _timing_flag and use public method that returns its value
  static int timing_flag;

  static Bool_t time_numInts();
  static void set_time_numInts(Bool_t flag);

private:
  static Bool_t _time_numInts;

protected:

  static RooTrace* _instance ;

  void dump3(std::ostream&, Bool_t sinceMarked) ;

  void mark3() ;
  void printObjectCounts3() ;
   
  void create2(const TObject* obj) ;
  void destroy2(const TObject* obj) ;

  void create3(const TObject* obj) ;
  void destroy3(const TObject* obj) ;

  void createSpecial3(const char* name, int size) ;
  void destroySpecial3(const char* name) ;

  void addPad(const TObject* ref, Bool_t doPad) ;
  Bool_t removePad(const TObject* ref) ;

  Bool_t _active ;
  Bool_t _verbose ;
  RooLinkedList _list ;
  RooLinkedList _markList ;
  std::map<TClass*,int> _objectCount ;
  std::map<std::string,int> _specialCount ;
  std::map<std::string,int> _specialSize ;

  ClassDef(RooTrace,1) // Memory tracer utility for RooFit objects
};


#endif
