// @(#)root/proofplayer:$Id$
// Author: G.Ganis July 2011

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


/** \class TProofMonSenderSQL
\ingroup proofkernel

TProofMonSender implementation for the SQL writers

*/

#include "TProofMonSenderSQL.h"

#include "TDSet.h"
#include "TFileInfo.h"
#include "THashList.h"
#include "TList.h"
#include "TUrl.h"
#include "TPluginManager.h"
#include "TProofDebug.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TObjString.h"
#include "TVirtualMonitoring.h"

////////////////////////////////////////////////////////////////////////////////
/// Main constructor

TProofMonSenderSQL::TProofMonSenderSQL(const char *serv, const char *user,
                                       const char *pass, const char *table,
                                       const char *dstab, const char *filestab)
                  : TProofMonSender(serv,"ProofMonSenderSQL"),
                    fDSetSendOpts("bulk,table=proofquerydsets"),
                    fFilesSendOpts("bulk,table=proofqueryfiles")
{
   fWriter = 0;
   // Init the sender instance using the plugin manager
   TPluginHandler *h = 0;
   if ((h = gROOT->GetPluginManager()->FindHandler("TVirtualMonitoringWriter", "SQL"))) {
      if (h->LoadPlugin() != -1) {
         fWriter = (TVirtualMonitoringWriter *) h->ExecPlugin(4, serv, user, pass, table);
         if (fWriter && fWriter->IsZombie()) SafeDelete(fWriter);
      }
   }
   // Flag this instance as valid if the writer initialization succeeded
   if (fWriter) ResetBit(TObject::kInvalidObject);

   // Set default send control options
   SetBit(TProofMonSender::kSendSummary);
   SetBit(TProofMonSender::kSendDataSetInfo);
   SetBit(TProofMonSender::kSendFileInfo);
   fSummaryVrs = 2;
   fDataSetInfoVrs = 1;
   fFileInfoVrs = 1;

   // Transfer verbosity requirements
   PDB(kMonitoring,1) if (fWriter) fWriter->Verbose(kTRUE);

   // Reformat the send options strings, if needed
   if (dstab && strlen(dstab) > 0) fDSetSendOpts.Form("bulk,table=%s", dstab);
   if (filestab && strlen(filestab) > 0) fFilesSendOpts.Form("bulk,table=%s", filestab);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TProofMonSenderSQL::~TProofMonSenderSQL()
{
   SafeDelete(fWriter);
}

////////////////////////////////////////////////////////////////////////////////
/// Send 'summary' record for the table 'proofquerylog'.
///
/// There are three versions of this record, corresponding the evolution
/// in time of the monitoring requirements.
///
/// The default version 2 corresponds to the table created with the following command:
///
/// CREATE TABLE proofquerylog (
///    id int(11) NOT NULL auto_increment,
///    proofuser varchar(32) NOT NULL,
///    proofgroup varchar(32) default NULL,
///    querybegin datetime default NULL,
///    queryend datetime default NULL,
///    walltime int(11) default NULL,
///    cputime float default NULL,
///    bytesread bigint(20) default NULL,
///    events bigint(20) default NULL,
///    totevents bigint(20) default NULL,
///    workers int(11) default NULL,
///    querytag varchar(64) NOT NULL,
///    vmemmxw bigint(20) default NULL,
///    rmemmxw bigint(20) default NULL,
///    vmemmxm bigint(20) default NULL,
///    rmemmxm bigint(20) default NULL,
///    numfiles int(11) default NULL,
///    missfiles int(11) default NULL,
///    status int(11) default NULL,
///    rootver varchar(32) NOT NULL,
///    PRIMARY KEY (id) );
///
/// Version 1 corresponds to the table created with the following command:
///    ('user','begin','end' instead of 'proofuser', 'querybegin', 'queryend';
///     no 'status', 'missfiles', 'rootver'; 'dataset' field with name(s) of
///     processed dataset(s))
///
/// CREATE TABLE proofquerylog (
///    id int(11) NOT NULL auto_increment,
///    user varchar(32) NOT NULL,
///    proofgroup varchar(32) default NULL,
///    begin datetime default NULL,
///    end datetime default NULL,
///    walltime int(11) default NULL,
///    cputime float default NULL,
///    bytesread bigint(20) default NULL,
///    events bigint(20) default NULL,
///    totevents bigint(20) default NULL,
///    workers int(11) default NULL,
///    querytag varchar(64) NOT NULL,
///    vmemmxw bigint(20) default NULL,
///    rmemmxw bigint(20) default NULL,
///    vmemmxm bigint(20) default NULL,
///    rmemmxm bigint(20) default NULL,
///    numfiles int(11) default NULL,
///    dataset varchar(512) NOT NULL,
///    PRIMARY KEY (id) );

Int_t TProofMonSenderSQL::SendSummary(TList *recs, const char *dumid)
{
   //
   // Version 0 corresponds to the table created with the following command:
   //    ('group' instead of 'proofgroup'; no 'querytag', 'vmemmxw',
   //     'rmemmxw', 'vmemmxm', 'rmemmxm', 'numfiles', 'dataset')
   //
   // CREATE TABLE proofquerylog (
   //    id int(11) NOT NULL auto_increment,
   //    user varchar(32) NOT NULL,
   //    group varchar(32) default NULL,
   //    begin datetime default NULL,
   //    end datetime default NULL,
   //    walltime int(11) default NULL,
   //    cputime float default NULL,
   //    bytesread bigint(20) default NULL,
   //    events bigint(20) default NULL,
   //    totevents bigint(20) default NULL,
   //    workers int(11) default NULL,
   //    PRIMARY KEY (id) );
   //
   //  Return 0 on success, -1 on any failure.

   if (!IsValid()) {
      Error("SendSummary", "invalid instance: do nothing!");
      return -1;
   }

   // Are we requested to send this info?
   if (!TestBit(TProofMonSender::kSendSummary)) return 0;

   PDB(kMonitoring,1) Info("SendSummary", "preparing (qid: '%s')", dumid);

   // Make sure we have something to send
   if (!recs || (recs && recs->GetSize() <= 0)) {
      Error("SendSummary", "records list undefined or empty!");
      return -1;
   }
   TList *xrecs = recs;

   TObject *dsn = 0;
   TNamed *nm = 0;
   // We may need to correct some variable names first
   if (fSummaryVrs > 1) {
      if ((nm = (TNamed *) recs->FindObject("user"))) nm->SetName("proofuser");
      if ((nm = (TNamed *) recs->FindObject("begin"))) nm->SetName("querybegin");
      if ((nm = (TNamed *) recs->FindObject("end"))) nm->SetName("queryend");
      if ((dsn = recs->FindObject("dataset"))) recs->Remove(dsn);
   } else if (fSummaryVrs == 0) {
      // Only the first records
      xrecs = new TList;
      xrecs->SetOwner(kFALSE);
      TIter nxr(recs);
      TObject *o = 0;
      while ((o = nxr())) {
         if (!strcmp(o->GetName(), "vmemmxw")) break;
         xrecs->Add(o);
      }
   }

   PDB(kMonitoring,1) Info("SendSummary", "sending (%d entries)", xrecs->GetSize());

   // Now we are ready to send
   Bool_t rc = fWriter->SendParameters(xrecs, dumid);

   // Restore the "dataset" entry in the list
   if (fSummaryVrs > 1 && dsn && xrecs == recs) {
      TObject *num = recs->FindObject("numfiles");
      if (num)
         recs->AddBefore(num, dsn);
      else
         recs->Add(dsn);
   }
   if (xrecs != recs) SafeDelete(xrecs);

   // Done
   return (rc ? 0 : -1);
}

////////////////////////////////////////////////////////////////////////////////
/// Post information about the processed dataset(s). The information is taken
/// from the TDSet object 'dset' and integrated with the missing files
/// information in the list 'missing'. The string 'qid' is the uninque
/// ID of the query; 'begin' the starting time.
///
/// The record is formatted for the table 'proofquerydsets'.
///
/// There are two versions of this record, with or without the starting time.
/// The starting time could be looked up from the summary record, if available.
///
/// The default version 1 corresponds to the table created with the following command:
///
/// CREATE TABLE proofquerydsets (
///    id int(11) NOT NULL auto_increment,
///    dsn varchar(512) NOT NULL,
///    querytag varchar(64) NOT NULL,
///    querybegin datetime default NULL,
///    numfiles int(11) default NULL,
///    missfiles int(11) default NULL,
///    PRIMARY KEY  (id),
///    KEY ix_querytag (querytag) );
///
/// Version 0 corresponds to the table created with the following command:
///    (no 'querybegin')
///
/// CREATE TABLE proofquerydsets (
///    id int(11) NOT NULL auto_increment,
///    dsn varchar(512) NOT NULL,
///    querytag varchar(64) NOT NULL,
///    numfiles int(11) default NULL,
///    missfiles int(11) default NULL,
///    PRIMARY KEY  (id),
///    KEY ix_querytag (querytag) );
///
/// The information is posted with a bulk insert.
///
/// Returns 0 on success, -1 on failure.

Int_t TProofMonSenderSQL::SendDataSetInfo(TDSet *dset, TList *missing,
                                          const char *begin, const char *qid)
{
   if (!IsValid()) {
      Error("SendDataSetInfo", "invalid instance: do nothing!");
      return -1;
   }

   // Are we requested to send this info?
   if (!TestBit(TProofMonSender::kSendDataSetInfo)) return 0;

   // The query id (tag) must be given
   if (!qid || (qid && strlen(qid) <= 0)) {
      Error("SendDataSetInfo", "query id (tag) undefined!");
      return -1;
   }
   // The dataset must be given
   if (!dset) {
      Error("SendDataSetInfo", "TDSet object undefined! (qid: '%s')", qid);
      return -1;
   }

   PDB(kMonitoring,1) Info("SendDataSetInfo", "preparing (qid: '%s')", qid);

   TList plets;
   // Extract the information and save it into the relevant multiplets
   TString dss(dset->GetName()), ds;
   Ssiz_t from = 0;
   while ((dss.Tokenize(ds, from , "[,| ]"))) {
      // Create a new TDSetPlet and add it to the list
      plets.Add(new TDSetPlet(ds.Data(), dset));
   }

   // Now try to count the files
   TDSetPlet *plet = 0;
   TIter nxpl(&plets);
   TObject *o = 0;
   TDSetElement *e = 0, *ee = 0;
   TDSet *dsete = 0;
   TIter nxe(dset->GetListOfElements());
   TString dse;
   while ((o = nxe())) {
      if ((e = dynamic_cast<TDSetElement *>(o))) {
         dse = e->GetDataSet();
         if (!dse.IsNull()) {
            nxpl.Reset();
            while ((plet = (TDSetPlet *) nxpl())) {
               if (dse == plet->GetName()) {
                  plet->fFiles += 1;
                  break;
               }
            }
         }
      } else if ((dsete = dynamic_cast<TDSet *>(o))) {
         PDB(kMonitoring,1)
            Info("SendDataSetInfo", "dset '%s' (%d files)",
                                    o->GetName(), dsete->GetListOfElements()->GetSize());
         TIter nxee(dsete->GetListOfElements());
         while ((ee = (TDSetElement *) nxee())) {
            dse = ee->GetDataSet();
            if (!dse.IsNull()) {
               nxpl.Reset();
               while ((plet = (TDSetPlet *) nxpl())) {
                  if (dse == plet->GetName()) {
                     plet->fFiles += 1;
                     plet->fDSet = dsete;
                     break;
                  }
               }
            }
         }
      } else {
         Warning("SendDataSetInfo", "ignoring unknown element type: '%s'", o->ClassName());
      }
   }

   // Now try to include the missing files info
   if (missing) {
      TFileInfo *fi = 0;
      TIter nxm(missing);
      TString dsfi, fn;
      while ((fi = (TFileInfo *) nxm())) {
         dsfi = fi->GetTitle();
         if (!dsfi.IsNull() && dsfi != "TFileInfo") {
            nxpl.Reset();
            while ((plet = (TDSetPlet *) nxpl())) {
               if (dsfi == plet->GetName()) {
                  fn = fi->GetCurrentUrl()->GetUrl();
                  if (plet->fDSet && plet->fDSet->GetListOfElements() &&
                      !(plet->fDSet->GetListOfElements()->FindObject(fn))) plet->fFiles += 1;
                  plet->fMissing += 1;
                  break;
               }
            }
         }
      }
   }

   // Now we can prepare the input for SendParameters
   TList values;
   TString ent("dsn,querytag,querybegin,numfiles,missfiles");
   if (fDataSetInfoVrs == 0) ent = "dsn,querytag,numfiles,missfiles";
   values.Add(new TObjString(ent.Data()));
   nxpl.Reset();
   while ((plet = (TDSetPlet *) nxpl())) {
      if (fDataSetInfoVrs == 0)
         ent.Form("'%s','%s',%d,%d", plet->GetName(), qid, plet->fFiles, plet->fMissing);
      else
         ent.Form("'%s','%s','%s',%d,%d", plet->GetName(), qid, begin, plet->fFiles, plet->fMissing);
      values.Add(new TObjString(ent.Data()));
   }

   PDB(kMonitoring,1)
      Info("SendDataSetInfo", "sending (%d entries)", values.GetSize());

   // Now we are ready to send
   Bool_t rc = fWriter->SendParameters(&values, fDSetSendOpts);

   // Done
   return (rc ? 0 : -1);
}

////////////////////////////////////////////////////////////////////////////////
/// Post information about the requested files. The information is taken
/// from the TDSet object 'dset' and integrated with the missing files
/// information in the list 'missing'. The string 'qid' is the unique
/// ID of the query; 'begin' the starting time.
///
/// The record is formatted for the table 'proofqueryfiles'.
///
/// There are two versions of this record, with or without the starting time.
/// The starting time could be looked up from the summary record, if available.
///
/// The default version 1 corresponds to the table created with the following command:
///
/// CREATE TABLE proofqueryfiles (
///    id int(11) NOT NULL auto_increment,
///    lfn varchar(255) NOT NULL,
///    path varchar(2048) NOT NULL,
///    querytag varchar(64) NOT NULL,
///    querybegin datetime default NULL,
///    status enum('Ok','Failed') NOT NULL default 'Ok',
///    PRIMARY KEY  (id),
///    KEY ix_querytag (querytag) );
///
/// Version 0 corresponds to the table created with the following command:
///    (no 'querybegin')
///
/// CREATE TABLE proofqueryfiles (
///    id int(11) NOT NULL auto_increment,
///    lfn varchar(255) NOT NULL,
///    path varchar(2048) NOT NULL,
///    querytag varchar(64) NOT NULL,
///    status enum('Ok','Failed') NOT NULL default 'Ok',
///    PRIMARY KEY  (id),
///    KEY ix_querytag (querytag) );
///
/// The information is posted with a bulk insert.
///
/// Returns 0 on success, -1 on failure.

Int_t TProofMonSenderSQL::SendFileInfo(TDSet *dset, TList *missing,
                                       const char *begin, const char *qid)
{
   if (!IsValid()) {
      Error("SendFileInfo", "invalid instance: do nothing!");
      return -1;
   }

   // Are we requested to send this info?
   if (!TestBit(TProofMonSender::kSendFileInfo)) return 0;

   // The query id (tag) must be given
   if (!qid || (qid && strlen(qid) <= 0)) {
      Error("SendFileInfo", "query id (tag) undefined!");
      return -1;
   }
   // The dataset must be given
   if (!dset) {
      Error("SendFileInfo", "TDSet object undefined! (qid: '%s')", qid);
      return -1;
   }

   PDB(kMonitoring,1) Info("SendFileInfo", "preparing (qid: '%s')", qid);
   THashList hmiss;
   if (missing) {
      TIter nxfm(missing);
      TFileInfo *fi = 0;
      while ((fi = (TFileInfo *)nxfm())) {
         hmiss.Add(new TObjString(fi->GetCurrentUrl()->GetUrl()));
      }
      PDB(kMonitoring,2) hmiss.Print();
   }

   TList values;
   TString ent("lfn,path,querytag,querybegin,status");
   if (fFileInfoVrs == 0)  ent = "lfn,path,querytag,status";
   values.Add(new TObjString(ent.Data()));

   // Create the file-plets
   TObject *o = 0;
   TDSetElement *e = 0, *ee = 0;
   TDSet *dsete = 0;
   TIter nxe(dset->GetListOfElements());
   TString fne, status;
   while ((o = nxe())) {
      if ((e = dynamic_cast<TDSetElement *>(o))) {
         fne = e->GetName();
         // Try to determine the status
         status = "Ok";
         if (hmiss.FindObject(fne)) status = "Failed";
         if (fFileInfoVrs == 0)
            ent.Form("'%s','%s','%s','%s'", gSystem->BaseName(fne), gSystem->GetDirName(fne).Data(),
                     qid, status.Data());
         else
            ent.Form("'%s','%s','%s','%s','%s'", gSystem->BaseName(fne), gSystem->GetDirName(fne).Data(),
                     qid, begin, status.Data());
         values.Add(new TObjString(ent.Data()));
      } else if ((dsete = dynamic_cast<TDSet *>(o))) {
         PDB(kMonitoring,1)
            Info("SendFileInfo", "dset '%s' (%d files)",
                                 o->GetName(), dsete->GetListOfElements()->GetSize());
         TIter nxee(dsete->GetListOfElements());
         while ((ee = (TDSetElement *) nxee())) {
            fne = ee->GetName();
            // Try to determine the status
            status = "Ok";
            if (hmiss.FindObject(fne)) status = "Failed";
            if (fFileInfoVrs == 0)
               ent.Form("'%s','%s','%s','%s'", gSystem->BaseName(fne), gSystem->GetDirName(fne).Data(),
                        qid, status.Data());
            else
               ent.Form("'%s','%s','%s','%s','%s'", gSystem->BaseName(fne), gSystem->GetDirName(fne).Data(),
                        qid, begin, status.Data());
            values.Add(new TObjString(ent.Data()));
         }
      } else {
         Warning("SendFileInfo", "ignoring unknown element type: '%s'", o->ClassName());
      }
   }

   PDB(kMonitoring,1) Info("SendFileInfo", "sending (%d entries)", values.GetSize());

   // Now we are ready to send
   Bool_t rc = fWriter->SendParameters(&values, fFilesSendOpts);

   // Done
   return (rc ? 0 : -1);
}
