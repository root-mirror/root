/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooMsgService.h,v 1.2 2007/07/13 21:50:24 wouter Exp $
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
#ifndef ROO_MSG_SERVICE
#define ROO_MSG_SERVICE

#include "TObject.h"
#include <string>
#include <vector>
#include <stack>
#include <map>
#include "RooCmdArg.h"
#include "RooGlobalFunc.h"
class RooAbsArg ;
class RooWorkspace ;

// Shortcut definitions 
#define coutI(a) RooMsgService::instance().log(this,RooFit::INFO,RooFit::a) 
#define coutP(a) RooMsgService::instance().log(this,RooFit::PROGRESS,RooFit::a) 
#define coutW(a) RooMsgService::instance().log(this,RooFit::WARNING,RooFit::a) 
#define coutE(a) RooMsgService::instance().log(this,RooFit::ERROR,RooFit::a) 
#define coutF(a) RooMsgService::instance().log(this,RooFit::FATAL,RooFit::a) 

// Skip the message prefix
#define ccoutD(a) RooMsgService::instance().log(this,RooFit::DEBUG,RooFit::a,kTRUE) 
#define ccoutI(a) RooMsgService::instance().log(this,RooFit::INFO,RooFit::a,kTRUE) 
#define ccoutP(a) RooMsgService::instance().log(this,RooFit::PROGRESS,RooFit::a,kTRUE) 
#define ccoutW(a) RooMsgService::instance().log(this,RooFit::WARNING,RooFit::a,kTRUE) 
#define ccoutE(a) RooMsgService::instance().log(this,RooFit::ERROR,RooFit::a,kTRUE) 
#define ccoutF(a) RooMsgService::instance().log(this,RooFit::FATAL,RooFit::a,kTRUE) 

// Message from given object instead of "this"
#define oocoutI(o,a) RooMsgService::instance().log(o,RooFit::INFO,RooFit::a) 
#define oocoutP(o,a) RooMsgService::instance().log(o,RooFit::PROGRESS,RooFit::a) 
#define oocoutW(o,a) RooMsgService::instance().log(o,RooFit::WARNING,RooFit::a) 
#define oocoutE(o,a) RooMsgService::instance().log(o,RooFit::ERROR,RooFit::a) 
#define oocoutF(o,a) RooMsgService::instance().log(o,RooFit::FATAL,RooFit::a) 

// Message from given object instead of "this" and skip message prefix
#define ooccoutD(o,a) RooMsgService::instance().log(o,RooFit::DEBUG,RooFit::a,kTRUE) 
#define ooccoutI(o,a) RooMsgService::instance().log(o,RooFit::INFO,RooFit::a,kTRUE) 
#define ooccoutP(o,a) RooMsgService::instance().log(o,RooFit::PROGRESS,RooFit::a,kTRUE) 
#define ooccoutW(o,a) RooMsgService::instance().log(o,RooFit::WARNING,RooFit::a,kTRUE) 
#define ooccoutE(o,a) RooMsgService::instance().log(o,RooFit::ERROR,RooFit::a,kTRUE) 
#define ooccoutF(o,a) RooMsgService::instance().log(o,RooFit::FATAL,RooFit::a,kTRUE) 

#ifndef _WIN32
#define ANYDEBUG (RooMsgService::_debugCount>0)
#else
#define ANYDEBUG (RooMsgService::anyDebug())
#endif

#define dologD(a) (ANYDEBUG && RooMsgService::instance().isActive(this,RooFit::a,RooFit::DEBUG))
#define dologI(a) (RooMsgService::instance().isActive(this,RooFit::a,RooFit::INFO))
#define dologP(a) (RooMsgService::instance().isActive(this,RooFit::a,RooFit::PROGRESS))
#define dologW(a) (RooMsgService::instance().isActive(this,RooFit::a,RooFit::WARNING))
#define dologE(a) (RooMsgService::instance().isActive(this,RooFit::a,RooFit::ERROR))
#define dologF(a) (RooMsgService::instance().isActive(this,RooFit::a,RooFit::FATAL))

#define oodologD(o,a) (ANYDEBUG && RooMsgService::instance().isActive(o,RooFit::a,RooFit::DEBUG))
#define oodologI(o,a) (RooMsgService::instance().isActive(o,RooFit::a,RooFit::INFO))
#define oodologP(o,a) (RooMsgService::instance().isActive(o,RooFit::a,RooFit::PROGRESS))
#define oodologW(o,a) (RooMsgService::instance().isActive(o,RooFit::a,RooFit::WARNING))
#define oodologE(o,a) (RooMsgService::instance().isActive(o,RooFit::a,RooFit::ERROR))
#define oodologF(o,a) (RooMsgService::instance().isActive(o,RooFit::a,RooFit::FATAL))

// Shortcuts definitions with conditional execution of print expression -- USE WITH CAUTION 

#define cxcoutD(a) if (ANYDEBUG && RooMsgService::instance().isActive(this,RooFit::a,RooFit::DEBUG)) RooMsgService::instance().log(this,RooFit::DEBUG,RooFit::a) 
#define ccxcoutD(a) if (ANYDEBUG && RooMsgService::instance().isActive(this,RooFit::a,RooFit::DEBUG)) RooMsgService::instance().log(this,RooFit::DEBUG,RooFit::a,kTRUE) 
#define oocxcoutD(o,a) if (ANYDEBUG && RooMsgService::instance().isActive(o,RooFit::a,RooFit::DEBUG)) RooMsgService::instance().log(o,RooFit::DEBUG,RooFit::a) 
#define ooccxcoutD(o,a) if (ANYDEBUG && RooMsgService::instance().isActive(o,RooFit::a,RooFit::DEBUG)) RooMsgService::instance().log(o,RooFit::DEBUG,RooFit::a,kTRUE) 
#define cxcoutI(a) if (RooMsgService::instance().isActive(this,RooFit::a,RooFit::INFO)) RooMsgService::instance().log(this,RooFit::INFO,RooFit::a) 
#define ccxcoutI(a) if (RooMsgService::instance().isActive(this,RooFit::a,RooFit::INFO)) RooMsgService::instance().log(this,RooFit::INFO,RooFit::a,kTRUE) 
#define oocxcoutI(o,a) if (RooMsgService::instance().isActive(o,RooFit::a,RooFit::INFO)) RooMsgService::instance().log(o,RooFit::INFO,RooFit::a) 
#define ooccxcoutI(o,a) if (RooMsgService::instance().isActive(o,RooFit::a,RooFit::INFO)) RooMsgService::instance().log(o,RooFit::INFO,RooFit::a,kTRUE) 
#define cxcoutP(a) if (RooMsgService::instance().isActive(this,RooFit::a,RooFit::PROGRESS)) RooMsgService::instance().log(this,RooFit::PROGRESS,RooFit::a)
#define ccxcoutP(a) if (RooMsgService::instance().isActive(this,RooFit::a,RooFit::PROGRESS)) RooMsgService::instance().log(this,RooFit::PROGRESS,RooFit::a,kTRUE)
#define oocxcoutP(o,a) if (RooMsgService::instance().isActive(o,RooFit::a,RooFit::PROGRESS)) RooMsgService::instance().log(o,RooFit::PROGRESS,RooFit::a)
#define ooccxcoutP(o,a) if (RooMsgService::instance().isActive(o,RooFit::a,RooFit::PROGRESS)) RooMsgService::instance().log(o,RooFit::PROGRESS,RooFit::a,kTRUE)
#define cxcoutW(a) if (RooMsgService::instance().isActive(this,RooFit::a,RooFit::WARNING)) RooMsgService::instance().log(this,RooFit::WARNING,RooFit::a) 
#define ccxcoutW(a) if (RooMsgService::instance().isActive(this,RooFit::a,RooFit::WARNING)) RooMsgService::instance().log(this,RooFit::WARNING,RooFit::a,kTRUE) 
#define oocxcoutW(o,a) if (RooMsgService::instance().isActive(o,RooFit::a,RooFit::WARNING)) RooMsgService::instance().log(o,RooFit::WARNING,RooFit::a) 
#define ooccxcoutW(o,a) if (RooMsgService::instance().isActive(o,RooFit::a,RooFit::WARNING)) RooMsgService::instance().log(o,RooFit::WARNING,RooFit::a,kTRUE) 
#define cxcoutE(a) if (RooMsgService::instance().isActive(this,RooFit::a,RooFit::ERROR)) RooMsgService::instance().log(this,RooFit::ERROR,RooFit::a) 
#define ccxcoutE(a) if (RooMsgService::instance().isActive(this,RooFit::a,RooFit::ERROR)) RooMsgService::instance().log(this,RooFit::ERROR,RooFit::a,kTRUE) 
#define oocxcoutE(o,a) if (RooMsgService::instance().isActive(o,RooFit::a,RooFit::ERROR)) RooMsgService::instance().log(o,RooFit::ERROR,RooFit::a)
#define ooccxcoutE(o,a) if (RooMsgService::instance().isActive(o,RooFit::a,RooFit::ERROR)) RooMsgService::instance().log(o,RooFit::ERROR,RooFit::a,kTRUE) 
#define cxcoutF(a) if (RooMsgService::instance().isActive(this,RooFit::a,RooFit::FATAL)) RooMsgService::instance().log(this,RooFit::FATAL,RooFit::a) 
#define ccxcoutF(a) if (RooMsgService::instance().isActive(this,RooFit::a,RooFit::FATAL)) RooMsgService::instance().log(this,RooFit::FATAL,RooFit::a,kTRUE) 
#define oocxcoutF(o,a) if (RooMsgService::instance().isActive(o,RooFit::a,RooFit::FATAL)) RooMsgService::instance().log(o,RooFit::FATAL,RooFit::a) 
#define ooccxcoutF(o,a) if (RooMsgService::instance().isActive(o,RooFit::a,RooFit::FATAL)) RooMsgService::instance().log(o,RooFit::FATAL,RooFit::a,kTRUE)

class RooMsgService : public TObject {
public:

  virtual ~RooMsgService() ;

  struct StreamConfig {
    public:

    void addTopic(RooFit::MsgTopic newTopic) {
      topic |= newTopic ;
    }

    void removeTopic(RooFit::MsgTopic oldTopic) {
      topic &= ~oldTopic ;
    }


    friend class RooMsgService ;
    
    Bool_t match(RooFit::MsgLevel level, RooFit::MsgTopic facility, const RooAbsArg* obj) ;
    Bool_t match(RooFit::MsgLevel level, RooFit::MsgTopic facility, const TObject* obj) ;
    
    Bool_t active ;
    Bool_t universal ;
    
    RooFit::MsgLevel minLevel ;
    Int_t    topic ;
    std::string objectName ;
    std::string className ;
    std::string baseClassName ;
    std::string tagName ;
    Color_t color ;
    Bool_t prefix ;
    
    std::ostream* os ;

  } ;

  // Access to instance
  static RooMsgService& instance();
  static Bool_t anyDebug() ;

  // User interface -- Add or delete reporting streams ;
  Int_t addStream(RooFit::MsgLevel level, const RooCmdArg& arg1=RooCmdArg(), const RooCmdArg& arg2=RooCmdArg(), const RooCmdArg& arg3=RooCmdArg(),
                    		  const RooCmdArg& arg4=RooCmdArg(), const RooCmdArg& arg5=RooCmdArg(), const RooCmdArg& arg6=RooCmdArg()); 
  void deleteStream(Int_t id) ;
  StreamConfig& getStream(Int_t id) { return _streams[id] ; }

  Int_t numStreams() const { return _streams.size() ; }
  void setStreamStatus(Int_t id, Bool_t active) ;
  Bool_t getStreamStatus(Int_t id) const ;

  void reset();

  void setGlobalKillBelow(RooFit::MsgLevel level) { _globMinLevel = level ; }
  RooFit::MsgLevel globalKillBelow() const { return _globMinLevel ; }

  void Print(Option_t *options= 0) const ;
  void showPid(Bool_t flag) { _showPid = flag ; }

  // Back end -- Send message or check if particular logging configuration is active
  std::ostream& log(const RooAbsArg* self, RooFit::MsgLevel level, RooFit::MsgTopic facility, Bool_t forceSkipPrefix=kFALSE) ;
  std::ostream& log(const TObject* self, RooFit::MsgLevel level, RooFit::MsgTopic facility, Bool_t forceSkipPrefix=kFALSE) ;
  Bool_t isActive(const RooAbsArg* self, RooFit::MsgTopic facility, RooFit::MsgLevel level) ;
  Bool_t isActive(const TObject* self, RooFit::MsgTopic facility, RooFit::MsgLevel level) ;

  static Int_t _debugCount ;
  std::map<int,std::string> _levelNames ;
  std::map<int,std::string> _topicNames ;

  // Print level support for RooFit-related messages that are not routed through RooMsgService (such as Minuit printouts)
  Bool_t silentMode() const { return _silentMode ; }  
  void setSilentMode(Bool_t flag) { _silentMode = flag ; }

  Int_t errorCount() const { return _errorCount ; }
  void clearErrorCount() { _errorCount = 0 ; }

  void saveState() ; 
  void restoreState() ;
  
  RooWorkspace* debugWorkspace() ;

  Int_t& debugCode() { return _debugCode ; }

protected:

  Int_t activeStream(const RooAbsArg* self, RooFit::MsgTopic facility, RooFit::MsgLevel level) ;
  Int_t activeStream(const TObject* self, RooFit::MsgTopic facility, RooFit::MsgLevel level) ;

  std::vector<StreamConfig> _streams ;
  std::stack<std::vector<StreamConfig> > _streamsSaved ;
  std::ostream* _devnull ;

  std::map<std::string,std::ostream*> _files ;
  RooFit::MsgLevel _globMinLevel ;
  RooFit::MsgLevel _lastMsgLevel ;

  Bool_t _silentMode ; 
  Bool_t _showPid ;

  Int_t _errorCount ;  

  // Private ctor -- singleton class
  RooMsgService() ;
  RooMsgService(const RooMsgService&) ;

  RooWorkspace* _debugWorkspace ;

  Int_t _debugCode ;
  
  ClassDef(RooMsgService,0) // RooFit Message Service Singleton class
};

#endif
