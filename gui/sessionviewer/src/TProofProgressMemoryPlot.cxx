// @(#)root/sessionviewer:$Id$
// Author: Anna Kreshuk 18/07/2008

/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TProofProgressMemoryPlot                                             //
//                                                                      //
// This class implements a dialog, used to display the memory footprint //
// on the proof workers and master. For the workers, memory is plotted  //
// as a function of number of events processed. For the master, it is   //
// plotted as a function of number of objects merged                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TProofProgressMemoryPlot.h"
#include "TProofProgressDialog.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TGListBox.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGMenu.h"
#include "TProofLog.h"
#include "TUrl.h"
#include "TProof.h"
#include "TError.h"
#include "TGFrame.h"
#include "TMacro.h"
#include "TObjString.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "snprintf.h"

#define kMemValuePos 8
#define kMemValuePosMaster 8
#define kEventNumberPos 13

ClassImp(TProofProgressMemoryPlot);

////////////////////////////////////////////////////////////////////////////////
/// Main constructor

TProofProgressMemoryPlot::TProofProgressMemoryPlot(TProofProgressDialog *d,
                                                   Int_t w, Int_t h)
                         : TGTransientFrame(gClient->GetRoot(),
                                            gClient->GetRoot(), w, h)
{
   fDialog = d;
   fProofLog = 0;
   fWPlot = 0;
   fMPlot = 0;
   fAPlot = 0;
   fFullLogs = kFALSE;

   // use hierarchical cleaning
   SetCleanup(kDeepCleanup);
   TGHorizontalFrame *htotal = new TGHorizontalFrame(this, w, h);
   //The frame for choosing workers
   TGVerticalFrame *vworkers = new TGVerticalFrame(htotal);
   TGLabel *label1 = new TGLabel(vworkers,"Choose workers:");

   //The list of workers
   fWorkers = BuildLogList(vworkers);
   fWorkers->Resize(102,52);
   fWorkers->SetMultipleSelections(kTRUE);

   //The SelectAll/ClearAll button
   TGPopupMenu *pm = new TGPopupMenu(gClient->GetRoot());
   pm->AddEntry("Select All", 0);
   pm->AddEntry("Clear All", 1);

   fAllWorkers = new TGSplitButton(vworkers, new TGHotString("Select ...            "), pm);
   fAllWorkers->Connect("ItemClicked(Int_t)", "TProofProgressMemoryPlot", this,
                        "Select(Int_t)");
   fAllWorkers->SetSplit(kFALSE);
   //select all for the first display
   Select(1);
   fWorkers->Select(0, kTRUE);
   fWorkers->Select(1, kTRUE);

   fPlot = new TGTextButton(vworkers, "Plot");
   fPlot->Connect("Clicked()", "TProofProgressMemoryPlot", this, "DoPlot()");
   vworkers->AddFrame(label1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 7, 2, 5, 2));
   vworkers->AddFrame(fAllWorkers, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 5, 2, 2, 2));
   vworkers->AddFrame(fWorkers, new TGLayoutHints(kLHintsExpandX | kLHintsTop | kLHintsExpandY, 2, 2, 5, 2));
   vworkers->AddFrame(fPlot, new TGLayoutHints(kLHintsExpandX | kLHintsBottom, 2, 2, 2, 2));
   htotal->AddFrame(vworkers, new TGLayoutHints(kLHintsCenterY | kLHintsLeft | kLHintsExpandY, 2, 2, 2, 2));

   fWorkersPlot = new TRootEmbeddedCanvas("WorkersCanvas", htotal, 300, 300);
   htotal->AddFrame(fWorkersPlot, new TGLayoutHints(kLHintsCenterY | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
   fMasterPlot = new TRootEmbeddedCanvas("MasterCanvas", htotal, 300, 300);
   htotal->AddFrame(fMasterPlot, new TGLayoutHints(kLHintsCenterY | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2 ,2));

   AddFrame(htotal, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
   TString title;
   title.Form("PROOF Memory Consumption: %s", (fDialog->fProof ?
              fDialog->fProof->GetMaster() : "<dummy>"));
   SetWindowName(title);
   SetIconName(title);

   MapSubwindows();
   Resize();
   CenterOnParent();
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TProofProgressMemoryPlot::~TProofProgressMemoryPlot()
{
   if (fProofLog){
      delete fProofLog;
      fProofLog = 0;
   }
   if (fMPlot){
      delete fMPlot;
      fMPlot = 0;
   }
   if (fWPlot){
      delete fWPlot;
      fWPlot = 0;
   }

   fProofLog = 0;
   fDialog->fMemWindow = 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Build the list of workers. For this, extract the logs and take the names
/// of TProofLogElements

TGListBox* TProofProgressMemoryPlot::BuildLogList(TGFrame *parent)
{
   TGListBox *c = new TGListBox(parent);
   c->AddEntry("average", 0);

   SafeDelete(fProofLog);
   fProofLog = 0;

   TProofMgr *mgr = TProof::Mgr(fDialog->fSessionUrl.Data());
   if (mgr) fProofLog = mgr->GetSessionLogs(0, 0, "Svc.*Memory");
   if (fDialog->fStatus==TProofProgressDialog::kRunning) {
      fFullLogs = kFALSE;
   } else {
      fFullLogs = kTRUE;
   }

   if (fProofLog) {
      TList *elem = fProofLog->GetListOfLogs();
      TIter next(elem);
      TProofLogElem *pe = 0;

      TString buf;
      Int_t is = 1;
      while ((pe=(TProofLogElem*)next())){
         TUrl url(pe->GetTitle());
         buf = TString::Format("%s %s", pe->GetName(), url.GetHost());
         c->AddEntry(buf.Data(), is);
         is++;
      }
   }
   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the canvases

void TProofProgressMemoryPlot::Clear(Option_t *)
{
   if (fWorkersPlot)
      fWorkersPlot->GetCanvas()->Clear();
   if (fMasterPlot)
      fMasterPlot->GetCanvas()->Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the plot from the logs

void TProofProgressMemoryPlot::DoPlot()
{
   Clear();

   if (!fProofLog || !fFullLogs ||
      (fDialog && fDialog->fStatus == TProofProgressDialog::kRunning)){

      SafeDelete(fProofLog);
      if (fDialog) {
         TProofMgr *mgr = TProof::Mgr(fDialog->fSessionUrl.Data());
         if (mgr) fProofLog = mgr->GetSessionLogs(0, 0, "Svc.*Memory");
         if (fDialog->fStatus==TProofProgressDialog::kRunning) {
            fFullLogs = kFALSE;
         } else {
            fFullLogs = kTRUE;
         }
      } else {
         Error("DoPlot", "TProofProgessDialog instance undefined - protocol error?");
         return;
      }
   }

   // Make sure we have something to parse
   if (!fProofLog) {
      Error("DoPlot", "could not get a valid instance of TProofLog");
      return;
   }

   char name[512]; //should be long enough

   TList *elem = fProofLog->GetListOfLogs();
   if (!elem) {Error("DoPlot", "No log elements\n"); return;}
   TIter next(elem);
   TProofLogElem *ple=0;

   Int_t iwelem = 0;
   Int_t imelem = 0;
   TGraph *gr=0;

   TList *selected = new TList;
   fWorkers->GetSelectedEntries(selected);
   TIter nextworker(selected);
   TGTextLBEntry *selworker;
   TLegend *legw = 0;
   TLegend *legm = 0;

   //delete the current multigraphs
   if (fWPlot){
      delete fWPlot;
      fWPlot = 0;
   }
   if (fMPlot) {
      delete fMPlot;
      fMPlot = 0;
   }

   //loop over the selected workers in the list
   Int_t max = -1;
   Int_t min = -1;
   while ((selworker=(TGTextLBEntry*)nextworker())){

      snprintf(name, sizeof(name)-1, "%s", selworker->GetText()->GetString());
      char *token;
      token = strtok(name, " ");
      if (token && !strcmp(token, "average")) { //change that to id comparison later
         gr = DoAveragePlot(max, min);
         if (gr && gr->GetN()>0){
            if (!fWPlot) {
               fWPlot = new TMultiGraph();
            }
            if (!legw) {
               legw = new TLegend(0.1, 0.7, 0.4, 0.9);
               legw->SetHeader("Workers");
            }
            gr->SetMarkerColor(1);
            gr->SetMarkerStyle(2);
            gr->SetMarkerSize(1);
            gr->SetLineWidth(2);
            gr->SetLineColor(1);
            fWPlot->Add(gr, "l");
            legw->AddEntry(gr, token, "l");
         }
         TProofLogElem *pltemp = (TProofLogElem*)elem->At(min+1);
         gr = DoWorkerPlot(pltemp);
         if (gr && gr->GetN()>0){
            if (!fWPlot) {
               fWPlot = new TMultiGraph();
            }
            if (!legw) {
               legw = new TLegend(0.1, 0.7, 0.4, 0.9);
               legw->SetHeader("Workers");
            }
            gr->SetLineWidth(2);
            gr->SetLineColor(2);
            gr->SetLineStyle(3);
            fWPlot->Add(gr, "l");
            legw->AddEntry(gr, TString::Format("%s - min", pltemp->GetName()) , "l");
         }
         pltemp = (TProofLogElem*)elem->At(max+1);
         gr = DoWorkerPlot(pltemp);
         if (gr && gr->GetN()>0){
            if (!fWPlot) {
               fWPlot = new TMultiGraph();
            }
            if (!legw) {
               legw = new TLegend(0.1, 0.7, 0.4, 0.9);
               legw->SetHeader("Workers");
            }
            gr->SetLineWidth(2);
            gr->SetLineColor(2);
            gr->SetLineStyle(2);
            fWPlot->Add(gr, "l");
            legw->AddEntry(gr, TString::Format("%s - max", pltemp->GetName()), "l");
         }

         continue;
      }


      ple = (TProofLogElem*)elem->FindObject(token);
      const char *role = ple->GetRole();
      if (role[0]=='w'){
         //role should be equal to "worker", only check the 1st char

         gr = DoWorkerPlot(ple);
         if (gr && gr->GetN()>0) {
            if (!fWPlot) {
               fWPlot = new TMultiGraph();
            }
            if (!legw) {
               legw = new TLegend(0.1, 0.7, 0.4, 0.9);
               legw->SetHeader("Workers");
            }
            gr->SetLineWidth(2);
            gr->SetLineColor(iwelem+3);
            fWPlot->Add(gr, "l");
            legw->AddEntry(gr, token, "l");
            iwelem++;
         }
      } else {
         //a master or submaster log
         //display without meaningful labels for now
         gr = DoMasterPlot(ple);
         if (gr && gr->GetN()>0){
            if (!fMPlot){
               fMPlot = new TMultiGraph();
            }
            if (!legm) {
               legm = new TLegend(0.1, 0.7, 0.4, 0.9);
               legm->SetHeader("Master");
            }
            gr->SetLineWidth(2);
            gr->SetLineColor(imelem+1);
            fMPlot->Add(gr, "l");
            legm->AddEntry(gr, token, "l");
            imelem++;
         }
      }

   }

   if (fWPlot){
      fWorkersPlot->GetCanvas()->cd();
      fWPlot->Draw("a");
      if (fWPlot->GetXaxis())
         fWPlot->GetXaxis()->SetTitle("Events Processed");
      if (fWPlot->GetYaxis())
         fWPlot->GetYaxis()->SetTitle("MBytes");
      if (legw) legw->Draw();

   }
   if (fMPlot) {
      fMasterPlot->GetCanvas()->cd();
      fMPlot->Draw("a");
      if (fMPlot->GetXaxis())
         fMPlot->GetXaxis()->SetTitle("Objects Merged");
      if (fMPlot->GetYaxis())
         fMPlot->GetYaxis()->SetTitle("MBytes");
      if (legm) legm->Draw();
   }
   fWorkersPlot->GetCanvas()->Update();
   fMasterPlot->GetCanvas()->Update();

   delete selected;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the average plots

TGraph *TProofProgressMemoryPlot::DoAveragePlot(Int_t &max_el, Int_t &min_el)
{
   TList *elem = fProofLog->GetListOfLogs();
   if (!elem) {
      Error("DoAveragePlot", "Empty log");
      return 0;
   }
   TIter next(elem);

   TProofLogElem *ple=0;
   Double_t max_av = 0;
   Double_t min_av = 10E9;

   Long64_t maxevent = 0;
   Long64_t step = -1;
   TObjString *curline = 0;
   TObjString *prevline = 0;
   Long64_t curevent_value;
   Long64_t prevevent_value;
   Long64_t *last = new Long64_t[elem->GetEntries()];
   Long64_t vmem = -1, rmem = -1, nevt = -1;
   TString token;
   Int_t ielem=0;
   for (Int_t i=0; i<elem->GetEntries(); i++) {
      last[i] = 0;
   }
   while ((ple = (TProofLogElem *)next())){
      //find the maximal entry processed in the last query
      const char *role = ple->GetRole();
      if (role[0] != 'w') continue; //skip the master log
      TList *lines = ple->GetMacro()->GetListOfLines();
      if (!lines || lines->GetSize() <= 0) continue;
      curline = (TObjString *) lines->Last();
      if (!curline) continue;
      curevent_value = 0;
      if (ParseLine(curline->String(), vmem, rmem, curevent_value) != 0) {
         Warning("DoAveragePlot", "error parsing line: '%s'", curline->String().Data());
         continue;
      }
      if (maxevent < curevent_value) maxevent = curevent_value;
      last[ielem] = curevent_value;
      if (step < 0) {
         // Find the step
         prevline = (TObjString *)lines->Before(curline);
         if (prevline) {
            prevevent_value = 0;
            if (ParseLine(prevline->String(), vmem, rmem, prevevent_value) != 0) {
               Warning("DoAveragePlot", "error parsing line: '%s'", curline->String().Data());
            } else {
               step = curevent_value - prevevent_value;
            }
         }
      }
      ielem++;
   }
   Int_t maxlines = Int_t(maxevent/(1.*step));
   //transform the array of last event numbers to an array of numbers of lines
   for (Int_t i=0; i<ielem; i++){
      last[i] /= step;
   }

   Double_t *av_mem = new Double_t[maxlines];
   Int_t *nw = new Int_t[maxlines];
   for (Int_t i=0; i<maxlines; i++){
      av_mem[i]=0;
      nw[i]=0;
   }
   next.Reset();
   ielem=0;
   Int_t iline=0;
   Double_t cur_av;
   while ((ple = (TProofLogElem*)next())){
      const char *role = ple->GetRole();
      if (role[0]!='w') continue;
      TList *lines = ple->GetMacro()->GetListOfLines();
      if (!lines || lines->GetSize() <= 0) continue;
      TIter prev(lines, kIterBackward);
      iline=0;
      cur_av = 0;
      while ((curline = (TObjString*)prev()) && iline<last[ielem]){
         // a backward loop, so that only the last query is counted
         vmem = 0;
         if (ParseLine(curline->String(), vmem, rmem, nevt) != 0) {
            Warning("DoWorkerPlot", "error parsing line: '%s'", curline->String().Data());
            continue;
         }
         av_mem[last[ielem] -1 - iline] += vmem; //last[ielem] is the number of lines for
         nw[last[ielem] -1 - iline]++;              //this query and this element
         if (last[ielem] > 0) cur_av += (Double_t)vmem / last[ielem];
         iline++;
      }
      if (cur_av > max_av){
         max_av = cur_av;
         max_el = ielem;
      }
      if (cur_av < min_av){
         min_av = cur_av;
         min_el = ielem;
      }
      ielem++;
   }

   TGraph *gr = new TGraph(maxlines);
   //make an average
   for (Int_t i=0; i<maxlines; i++){
      gr->SetPoint(i, (i+1)*step, av_mem[i]/(nw[i]*1024.));
   }
   delete [] av_mem;
   av_mem = 0;
   delete [] nw;
   nw = 0;
   delete [] last;
   last = 0;
   return gr;

}

////////////////////////////////////////////////////////////////////////////////
/// Extract from line 'l' the virtual memory 'v', the resident memory 'r' and the
/// number of events 'e'.
/// The line is assumed to be in the form
/// "... Memory 130868 virtual 31540 ... event 5550"
/// The fields are only filled if >= 0 .
/// Return 0 on success, -1 if any of the values coudl not be filled (the output
/// fields are not touched in such a case).

Int_t TProofProgressMemoryPlot::ParseLine(TString l,
                                          Long64_t &v, Long64_t &r, Long64_t &e)
{
   // Something to parse is mandatory
   if (l.IsNull()) return -1;

   // At least one field needs to be filled
   if (v < 0 && r < 0 && e < 0) return 0;

   // Position at the start of the relevant info
   Int_t from = kNPOS;
   if ((from = l.Index("Memory")) == kNPOS) return -1;

   // Prepare extraction
   from += 7;
   TString tok;

   // The virtual memory
   if (v >= 0) {
      if (!l.Tokenize(tok, from, " ")) return -1;
      v = tok.Atoll();
   }

   // The resident memory
   if (r >= 0) {
      if ((from = l.Index("virtual", from)) == kNPOS) return -1;
      from += 8;
      if (!l.Tokenize(tok, from, " ")) return -1;
      r = tok.Atoll();
   }

   // The number of events
   if (e >= 0) {
      if ((from = l.Index("event", from)) == kNPOS) return -1;
      from += 6;
      if (!l.Tokenize(tok, from, " ")) return -1;
      e = tok.Atoll();
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a memory consumption graph for a worker

TGraph *TProofProgressMemoryPlot::DoWorkerPlot(TProofLogElem *ple)
{
   TObjString *curline;
   TList *lines = ple->GetMacro()->GetListOfLines();
   if (!lines) {
      //the log is empty so far
      return 0;
   }

   Long64_t vmem = -1, rmem = -1, nevt = -1;

   //find the last event value
   curline = (TObjString*)lines->Last();
   Long64_t lastevent_value = 0;
   if (ParseLine(curline->String(), vmem, rmem, lastevent_value) != 0) {
      Error("DoWorkerPlot", "error parsing line: '%s'", curline->String().Data());
      return 0;
   }

   //find the step
   TObjString *prevline = (TObjString*)lines->Before(curline);
   Long64_t prevevent_value = 0;
   if (prevline && ParseLine(prevline->String(), vmem, rmem, prevevent_value) != 0) {
      Error("DoWorkerPlot", "error parsing line: '%s'", prevline->String().Data());
      return 0;
   }
   Long64_t step = lastevent_value - prevevent_value;
   if (step <= 0) {
      Error("DoWorkerPlot", "null or negative step (%lld) - cannot continue", step);
      return 0;
   }

   Int_t nlines = lastevent_value/step;
   TGraph *gr = new TGraph(nlines);

   TIter prevl(lines, kIterBackward);
   Int_t iline = 0;
   TString token;
   while ((curline = (TObjString*)prevl()) && iline<nlines){
      //iterate backwards so that only lines for the last query are taken
      vmem = 0;
      if (ParseLine(curline->String(), vmem, rmem, nevt) != 0) {
         Warning("DoWorkerPlot", "error parsing line: '%s'", curline->String().Data());
         continue;
      }
      gr->SetPoint(nlines-1-iline, lastevent_value-iline*step, vmem/1024.);
      iline++;
   }

   return gr;
}

////////////////////////////////////////////////////////////////////////////////
/// a master or submaster log
/// display without meaningful labels for now

TGraph *TProofProgressMemoryPlot::DoMasterPlot(TProofLogElem *ple)
{
   TList *lines = ple->GetMacro()->GetListOfLines();
   TIter prevline(lines, kIterBackward);
   Int_t iline=0;
   TObjString *curline;
   //count the number of lines that belong to the last query
   while ((curline = (TObjString*)prevline())) {
      if (curline->String().Contains("Start")) break;
      iline++;
   }

   Long64_t vmem = -1, rmem = -1, nevt = -1;

   Int_t nlines = iline;
   TString token;
   TGraph *gr = new TGraph(nlines);
   prevline.Reset();
   iline = 0;
   while ((curline = (TObjString*)prevline()) && iline<nlines) {
      //iterate backwards so that only lines for the last query are taken
      vmem = 0;
      if (ParseLine(curline->String(), vmem, rmem, nevt) != 0) {
         Warning("DoWorkerPlot", "error parsing line: '%s'", curline->String().Data());
         continue;
      }
      gr->SetPoint(nlines-iline, nlines-iline, vmem/1024.);
      iline++;
   }
   return gr;
}

////////////////////////////////////////////////////////////////////////////////
///actions of select all/clear all button

void TProofProgressMemoryPlot::Select(Int_t id)
{
   Int_t nen = fWorkers->GetNumberOfEntries();
   Bool_t sel = id ? 0 : 1;

   for (Int_t ie=0; ie<nen; ie++) {
      fWorkers->Select(ie, sel);
   }
}
