// @(#)root/base:$Id$
// Author: Andreas-Joachim Peters   15/05/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualMonitoring
#define ROOT_TVirtualMonitoring

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualMonitoring                                                   //
//                                                                      //
// Provides the interface for Monitoring plugins.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

#ifdef R__LESS_INCLUDES
class TList;
class TMap;
#else
#include "TList.h"
#include "TMap.h"
#endif

class TFile;

class TVirtualMonitoringWriter : public TNamed {

private:

   TVirtualMonitoringWriter(const TVirtualMonitoringWriter&) = delete;
   TVirtualMonitoringWriter& operator=(const TVirtualMonitoringWriter&) = delete;

   Double_t fValue;  // double monitor value

protected:
   TList     *fTmpOpenPhases;       // To store open phases when there is not yet an object

public:
   TVirtualMonitoringWriter() : TNamed(), fValue(0), fTmpOpenPhases(0) { }
   TVirtualMonitoringWriter(const char *name, Double_t value)
     : TNamed(name, ""), fValue(value), fTmpOpenPhases(0) { }

   virtual ~TVirtualMonitoringWriter();

   // TFile related info. In general they are gathered and sent only sometimes as summaries
   virtual Bool_t SendFileCloseEvent(TFile * /*file*/)
      { MayNotUse("SendFileCloseEvent"); return kFALSE; }
   virtual Bool_t SendFileReadProgress(TFile * /*file*/)
      { MayNotUse("SendFileReadProgress"); return kFALSE; }
   virtual Bool_t SendFileWriteProgress(TFile * /*file*/)
      { MayNotUse("SendFileWriteProgress"); return kFALSE; }

   virtual Bool_t SendParameters(TList * /*valuelist*/, const char * /*identifier*/ = 0)
      { MayNotUse("SendParameters"); return kFALSE; }
   virtual Bool_t SendInfoTime() { MayNotUse("SendInfoTime"); return kFALSE; }
   virtual Bool_t SendInfoUser(const char * /*user*/ = 0) { MayNotUse("SendInfoUser"); return kFALSE; }
   virtual Bool_t SendInfoDescription(const char * /*jobtag*/) { MayNotUse("SendInfoDescription"); return kFALSE; }
   virtual Bool_t SendInfoStatus(const char * /*status*/) { MayNotUse("SendInfoStatus"); return kFALSE; }

   // An Open might have several phases, and the timings might be interesting
   // to report
   // The info is only gathered into openphasestime, and sent when forcesend=kTRUE
   virtual Bool_t SendFileOpenProgress(TFile * /*file*/, TList * /*openphases*/,
                                       const char * /*openphasename*/,
                                       Bool_t /*forcesend*/ = kFALSE )
      { MayNotUse("SendFileOpenProgress"); return kFALSE; }

   virtual Bool_t SendProcessingStatus(const char * /*status*/, Bool_t /*restarttimer*/ = kFALSE)
      { MayNotUse("SendProcessingStatus"); return kFALSE; }
   virtual Bool_t SendProcessingProgress(Double_t /*nevent*/, Double_t /*nbytes*/, Bool_t /*force*/ = kFALSE)
      { MayNotUse("SendProcessingProgress"); return kFALSE; }
   virtual void   SetLogLevel(const char * /*loglevel*/ = "WARNING")
      { MayNotUse("SetLogLevel"); };
   virtual void   Verbose(Bool_t /*onoff*/) { MayNotUse("Verbose"); }

   ClassDef(TVirtualMonitoringWriter,0)  // ABC for Sending Monitoring Information
};


class TVirtualMonitoringReader : public TNamed {

public:
   TVirtualMonitoringReader( const char * /*serviceurl*/ ="") { }
   virtual ~TVirtualMonitoringReader() { }

   virtual void   DumpResult() { MayNotUse("DumpResult"); }
   virtual void   GetValues(const char * /*farmName*/, const char * /*clusterName*/,
                            const char * /*nodeName*/, const char * /*paramName*/,
                            Long_t /*min*/, Long_t /*max*/, Bool_t /*debug*/ = kFALSE)
      { MayNotUse("GetValues"); }
   virtual void   GetLastValues(const char * /*farmName*/, const char * /*clusterName*/,
                                const char * /*nodeName*/, const char * /*paramName*/,
                                Bool_t /*debug*/ =kFALSE)
      { MayNotUse("GetLastValues"); }
   virtual void   ProxyValues(const char * /*farmName*/, const char * /*clusterName*/,
                              const char * /*nodeName*/, const char * /*paramName*/,
                              Long_t /*min*/, Long_t /*max*/, Long_t /*lifetime*/)
      { MayNotUse("ProxyValues"); }

   virtual TMap *GetMap() { MayNotUse("GetMap"); return 0; }
   virtual void DeleteMap(TMap * /*map*/) { MayNotUse("DeleteMap"); }

   ClassDef(TVirtualMonitoringReader, 1) // ABC for Reading Monitoring Information
};


R__EXTERN TVirtualMonitoringWriter *gMonitoringWriter;
R__EXTERN TVirtualMonitoringReader *gMonitoringReader;


#endif
