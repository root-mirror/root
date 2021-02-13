// @(#)root/proofplayer:$Id$
// Author: G.Ganis July 2011

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


/** \class TProofMonSenderML
\ingroup proofkernel

TProofMonSender implementation for the ML writer

*/

#include "TProofMonSenderML.h"

#include "TDSet.h"
#include "TFileInfo.h"
#include "THashList.h"
#include "TList.h"
#include "TParameter.h"
#include "TPluginManager.h"
#include "TProofDebug.h"
#include "TROOT.h"
#include "TUrl.h"
#include "TSystem.h"
#include "TObjString.h"
#include "TVirtualMonitoring.h"

////////////////////////////////////////////////////////////////////////////////
/// Main constructor

TProofMonSenderML::TProofMonSenderML(const char *serv, const char *tag,
                                     const char *id, const char *subid,
                                     const char *opt)
                  : TProofMonSender(serv, "ProofMonSenderML")
{
   fWriter = 0;
   // Init the sender instance using the plugin manager
   TPluginHandler *h = 0;
   if ((h = gROOT->GetPluginManager()->FindHandler("TVirtualMonitoringWriter", "MonaLisa"))) {
      if (h->LoadPlugin() != -1) {
         fWriter = (TVirtualMonitoringWriter *) h->ExecPlugin(5, serv, tag, id, subid, opt);
         if (fWriter && fWriter->IsZombie()) SafeDelete(fWriter);
      }
   }
   // Flag this instance as valid if the writer initialization succeeded
   if (fWriter) ResetBit(TObject::kInvalidObject);

   // Set default send control options
   SetBit(TProofMonSender::kSendSummary);
   ResetBit(TProofMonSender::kSendDataSetInfo);
   ResetBit(TProofMonSender::kSendFileInfo);
   fSummaryVrs = 1;
   fDataSetInfoVrs = 1;
   fFileInfoVrs = 1;

   // Transfer verbosity requirements
   PDB(kMonitoring,1) if (fWriter) fWriter->Verbose(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TProofMonSenderML::~TProofMonSenderML()
{
   SafeDelete(fWriter);
}

////////////////////////////////////////////////////////////////////////////////
/// Send summary record
///
/// There are three versions of this record, corresponding the evolution
/// in time of the monitoring requirements.
///
/// The default version 2 contains the following information
///
///    user                  XRD_STRING
///    proofgroup            XRD_STRING
///    begin                 XRD_STRING
///    end                   XRD_STRING
///    walltime              XRD_REAL64
///    cputime               XRD_REAL64
///    bytesread             XRD_REAL64
///    events                XRD_REAL64
///    totevents             XRD_REAL64
///    workers               XRD_REAL64
///    vmemmxw               XRD_REAL64
///    rmemmxw               XRD_REAL64
///    vmemmxm               XRD_REAL64
///    rmemmxm               XRD_REAL64
///    numfiles              XRD_REAL64
///    missfiles             XRD_REAL64
///    status                XRD_REAL64
///    rootver               XRD_STRING
///
/// Version 1 contains the following information
///    (no 'status', 'missfiles', 'rootver'; 'dataset' field with name(s) of
///     processed dataset(s))
///
///    user                  XRD_STRING
///    proofgroup            XRD_STRING
///    begin                 XRD_STRING
///    end                   XRD_STRING
///    walltime              XRD_REAL64
///    cputime               XRD_REAL64
///    bytesread             XRD_REAL64
///    events                XRD_REAL64
///    totevents             XRD_REAL64
///    workers               XRD_REAL64
///    vmemmxw               XRD_REAL64
///    rmemmxw               XRD_REAL64
///    vmemmxm               XRD_REAL64
///    rmemmxm               XRD_REAL64
///    numfiles              XRD_REAL64
///    dataset               XRD_STRING
///
/// Version 0 contains the following information
///    ('group' instead of 'proofgroup'; no 'vmemmxw',
///     'rmemmxw', 'vmemmxm', 'rmemmxm', 'numfiles', 'dataset')
///
///    user                  XRD_STRING
///    group                 XRD_STRING
///    begin                 XRD_STRING
///    end                   XRD_STRING
///    walltime              XRD_REAL64
///    cputime               XRD_REAL64
///    bytesread             XRD_REAL64
///    events                XRD_REAL64
///    totevents             XRD_REAL64
///    workers               XRD_REAL64
///
///  Return 0 on success, -1 on any failure.

Int_t TProofMonSenderML::SendSummary(TList *recs, const char *id)
{
   if (!IsValid()) {
      Error("SendSummary", "invalid instance: do nothing!");
      return -1;
   }

   // Are we requested to send this info?
   if (!TestBit(TProofMonSender::kSendSummary)) return 0;

   // Make sure we have something to send
   if (!recs || (recs && recs->GetSize() <= 0)) {
      Error("SendSummary", "records list undefined or empty!");
      return -1;
   }
   TList *xrecs = recs;

   PDB(kMonitoring,1) Info("SendSummary", "preparing (qid: '%s')", id);

   // Do not send duplicated information
   TObject *qtag = recs->FindObject("querytag");
   if (qtag) recs->Remove(qtag);

   TObject *dsn = 0;
   // We may need to correct some variable names first
   if (fSummaryVrs > 1) {
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

   PDB(kMonitoring,1)
      Info("SendSummary", "sending (%d entries)", xrecs->GetSize());

   // Now we are ready to send
   Bool_t rc = fWriter->SendParameters(xrecs, id);

   // Restore the "dataset" entry in the list
   if (fSummaryVrs > 1 && dsn && xrecs == recs) {
      TObject *num = recs->FindObject("numfiles");
      if (num)
         recs->AddBefore(num, dsn);
      else
         recs->Add(dsn);
   }
   // Restore the "querytag" entry in the list
   if (qtag) {
      TObject *wrks = recs->FindObject("workers");
      if (wrks)
         recs->AddAfter(wrks, qtag);
      else
         recs->Add(qtag);
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
/// The records sent by this call will appear with ids 'dataset_<dataset_name_hash>'
///
/// There are two versions of this record, with or without the starting time.
/// The starting time could be looked up from the summary record, if available.
///
/// The default version 1 contains the following information:
///
///    dsn              XRD_STRING
///    querytag         XRD_STRING
///    querybegin       XRD_STRING
///    numfiles         XRD_REAL64
///    missfiles        XRD_REAL64
///
/// Version 0 contains the following information:
///    (no 'querybegin')
///
///    dsn              XRD_STRING
///    querytag         XRD_STRING
///    numfiles         XRD_REAL64
///    missfiles        XRD_REAL64
///
/// The information is posted with a bulk insert.
///
/// Returns 0 on success, -1 on failure.

Int_t TProofMonSenderML::SendDataSetInfo(TDSet *dset, TList *missing,
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

   PDB(kMonitoring,1)
      Info("SendDataSetInfo", "preparing (qid: '%s')", qid);

   TList plets;
   // Extract the information and save it into the relevant multiplets
   TString dss(dset->GetName()), ds;
   Ssiz_t from = 0;
   while ((dss.Tokenize(ds, from , "[,| ]"))) {
      // Create a new TDSetPlet and add it to the list
      plets.Add(new TDSetPlet(ds.Data()));
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

   // Prepare objects to be sent
   TList values;
   TNamed *nm_dsn = new TNamed("dsn", "");
   values.Add(nm_dsn);
   TNamed *nm_querytag = new TNamed("querytag", qid);
   values.Add(nm_querytag);
   TNamed *nm_begin = 0;
   if (fDataSetInfoVrs > 0) {
      nm_begin = new TNamed("begin", begin);
      values.Add(nm_begin);
   }
   TParameter<Int_t> *pi_numfiles = new TParameter<Int_t>("numfiles", -1);
   values.Add(pi_numfiles);
   TParameter<Int_t> *pi_missfiles = new TParameter<Int_t>("missfiles", -1);
   values.Add(pi_missfiles);

   PDB(kMonitoring,1)
      Info("SendDataSetInfo", "sending (%d entries)", plets.GetSize());

   Bool_t rc = kTRUE;
   TString dsnh;
   nxpl.Reset();
   while ((plet = (TDSetPlet *) nxpl())) {
      nm_dsn->SetTitle(plet->GetName());
      pi_numfiles->SetVal(plet->fFiles);
      pi_missfiles->SetVal(plet->fMissing);
      dsnh.Form("dataset_%x", TString(plet->GetName()).Hash());
      if (!(rc = fWriter->SendParameters(&values, dsnh.Data()))) break;
   }

   // Done
   return (rc ? 0 : -1);
}

////////////////////////////////////////////////////////////////////////////////
/// Post information about the requested files. The information is taken
/// from the TDSet object 'dset' and integrated with the missing files
/// information in the list 'missing'. The string 'qid' is the unique
/// ID of the query; 'begin' the starting time.
///
/// The records sent by this call will appear with ids 'file_<file_name_hash>'
///
/// There are two versions of this record, with or without the starting time.
/// The starting time could be looked up from the summary record, if available.
///
/// The default version 1 contains the following information:
///
///    lfn              XRD_STRING
///    path             XRD_STRING
///    querytag         XRD_STRING
///    querybegin       XRD_STRING
///    status           XRD_REAL64
///
/// Version 0 contains the following information:
///    (no 'querybegin')
///
///    lfn              XRD_STRING
///    path             XRD_STRING
///    querytag         XRD_STRING
///    status           XRD_REAL64
///
/// The information is posted with a bulk insert.
///
/// Returns 0 on success, -1 on failure.

Int_t TProofMonSenderML::SendFileInfo(TDSet *dset, TList *missing,
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
      hmiss.Print();
   }

   // Prepare objects to be sent
   TList values;
   TNamed *nm_lnf = new TNamed("lnf", "");
   values.Add(nm_lnf);
   TNamed *nm_path = new TNamed("path", "");
   values.Add(nm_path);
   TNamed *nm_querytag = new TNamed("querytag", qid);
   values.Add(nm_querytag);
   TNamed *nm_begin = 0;
   if (fFileInfoVrs > 0) {
      nm_begin = new TNamed("begin", begin);
      values.Add(nm_begin);
   }
   TParameter<Int_t> *pi_status = new TParameter<Int_t>("status", -1);
   values.Add(pi_status);

   PDB(kMonitoring,1)
      Info("SendFileInfo", "sending (%d entries)",
                           dset->GetListOfElements()->GetSize());

   // Loop over files
   Bool_t rc = kTRUE;
   TObject *o = 0;
   TDSetElement *e = 0, *ee = 0;
   TDSet *dsete = 0;
   TIter nxe(dset->GetListOfElements());
   TString fne, fneh;
   Int_t status = -1;
   while ((o = nxe())) {
      if ((e = dynamic_cast<TDSetElement *>(o))) {
         fne = e->GetName();
         // Try to determine the status
         status = 1;
         if (hmiss.FindObject(fne)) status = 0;
         // Prepare the parameters list
         nm_lnf->SetTitle(gSystem->BaseName(fne));
         nm_path->SetTitle(gSystem->GetDirName(fne));
         pi_status->SetVal(status);
         fneh.Form("file_%x", TString(TUrl(fne.Data()).GetFile()).Hash());
         if (!(rc = fWriter->SendParameters(&values, fneh.Data()))) break;
      } else if ((dsete = dynamic_cast<TDSet *>(o))) {
         PDB(kMonitoring,1)
            Info("SendFileInfo", "dset '%s' (%d files)",
                                 o->GetName(), dsete->GetListOfElements()->GetSize());
         TIter nxee(dsete->GetListOfElements());
         while ((ee = (TDSetElement *) nxee())) {
            fne = ee->GetName();
            // Try to determine the status
            status = 1;
            if (hmiss.FindObject(fne)) status = 0;
            // Prepare the parameters list
            nm_lnf->SetTitle(gSystem->BaseName(fne));
            nm_path->SetTitle(gSystem->GetDirName(fne));
            pi_status->SetVal(status);
            fneh.Form("file_%x", TString(TUrl(fne.Data()).GetFile()).Hash());
            if (!(rc = fWriter->SendParameters(&values, fneh.Data()))) break;
         }
      } else {
         Warning("SendFileInfo", "ignoring unknown element type: '%s'", o->ClassName());
      }
   }

   // Done
   return (rc ? 0 : -1);
}
