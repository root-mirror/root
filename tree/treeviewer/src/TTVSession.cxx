// @(#)root/treeviewer:$Id$
//Author : Andrei Gheata   21/02/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include <iostream>
#include <fstream>
#include "TTVSession.h"
#include "TTreeViewer.h"
#include "TTVLVContainer.h"
#include "TClonesArray.h"
#include "TInterpreter.h"
#include "snprintf.h"


ClassImp(TTVRecord);

/** \class TTVRecord
I/O classes for TreeViewer session handling.
*/

////////////////////////////////////////////////////////////////////////////////
/// TTVRecord default constructor

TTVRecord::TTVRecord()
{
   fName = "";
   fScanRedirected = kFALSE;
   fCutEnabled     = kTRUE;
   fUserCode = "";
   fAutoexec = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute user-defined code

void TTVRecord::ExecuteUserCode()
{
   if (fUserCode.Length()) {
      char code[250];
      code[0] = 0;
      snprintf(code,250, "%s", fUserCode.Data());
      gInterpreter->ProcessLine(code);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Populate members from treeviewer tv

void TTVRecord::FormFrom(TTreeViewer *tv)
{
   if (!tv)  return;
   fX        = tv->ExpressionItem(0)->GetTrueName();
   fXAlias   = tv->ExpressionItem(0)->GetAlias();
   fY        = tv->ExpressionItem(1)->GetTrueName();
   fYAlias   = tv->ExpressionItem(1)->GetAlias();
   fZ        = tv->ExpressionItem(2)->GetTrueName();
   fZAlias   = tv->ExpressionItem(2)->GetAlias();
   fCut      = tv->ExpressionItem(3)->GetTrueName();
   fCutAlias = tv->ExpressionItem(3)->GetAlias();
   fOption   = tv->GetGrOpt();
   fScanRedirected = tv->IsScanRedirected();
   fCutEnabled = tv->IsCutEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Change treeviewer status to this record

void TTVRecord::PlugIn(TTreeViewer *tv)
{
   TTVLVEntry *item;
   // change X expression
   item = tv->ExpressionItem(0);
   item->SetExpression(fX.Data(), fXAlias.Data());
   item = tv->ExpressionItem(1);
   item->SetExpression(fY.Data(), fYAlias.Data());
   item = tv->ExpressionItem(2);
   item->SetExpression(fZ.Data(), fZAlias.Data());
   item = tv->ExpressionItem(3);
   item->SetExpression(fCut.Data(), fCutAlias.Data());
   tv->SetGrOpt(fOption.Data());
   tv->SetScanRedirect(fScanRedirected);
   tv->SetCutMode(fCutEnabled);
   if (fCutEnabled)
      item->SetSmallPic(gClient->GetPicture("cut_t.xpm"));
   else
      item->SetSmallPic(gClient->GetPicture("cut-disable_t.xpm"));
}

////////////////////////////////////////////////////////////////////////////////
/// Save the TTVRecord in a C++ macro file

void TTVRecord::SaveSource(std::ofstream &out)
{
   char quote = '"';
   out <<"//--- tree viewer record"<<std::endl;
   out <<"   tv_record = tv_session->AddRecord(kTRUE);"<<std::endl;
   out <<"   tv_session->SetRecordName("<<quote<<GetName()<<quote<<");"<<std::endl;
   out <<"   tv_record->fX        = "<<quote<<fX.Data()<<quote<<";"<<std::endl;
   out <<"   tv_record->fY        = "<<quote<<fY.Data()<<quote<<";"<<std::endl;
   out <<"   tv_record->fZ        = "<<quote<<fZ.Data()<<quote<<";"<<std::endl;
   out <<"   tv_record->fCut      = "<<quote<<fCut.Data()<<quote<<";"<<std::endl;
   out <<"   tv_record->fXAlias   = "<<quote<<fXAlias.Data()<<quote<<";"<<std::endl;
   out <<"   tv_record->fYAlias   = "<<quote<<fYAlias.Data()<<quote<<";"<<std::endl;
   out <<"   tv_record->fZAlias   = "<<quote<<fZAlias.Data()<<quote<<";"<<std::endl;
   out <<"   tv_record->fCutAlias = "<<quote<<fCutAlias.Data()<<quote<<";"<<std::endl;
   out <<"   tv_record->fOption   = "<<quote<<fOption.Data()<<quote<<";"<<std::endl;
   if (fScanRedirected)
      out <<"   tv_record->fScanRedirected = kTRUE;"<<std::endl;
   else
      out <<"   tv_record->fScanRedirected = kFALSE;"<<std::endl;
   if (fCutEnabled)
      out <<"   tv_record->fCutEnabled = kTRUE;"<<std::endl;
   else
      out <<"   tv_record->fCutEnabled = kFALSE;"<<std::endl;
   if (fUserCode.Length()) {
      out <<"   tv_record->SetUserCode(\""<<fUserCode.Data()<<"\");"<<std::endl;
      if (fAutoexec) {
         out <<"   tv_record->SetAutoexec();"<<std::endl;
      }
   }
}

ClassImp(TTVSession);

/** \class TTVSession
I/O classes for TreeViewer session handling.
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TTVSession::TTVSession(TTreeViewer *tv)
{
   fName    = "";
   fList    = new TClonesArray("TTVRecord", 100); // is 100 enough ?
   fViewer  = tv;
   fCurrent = 0;
   fRecords = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TTVSession::~TTVSession()
{
   fList->Delete();
   delete fList;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a record

TTVRecord *TTVSession::AddRecord(Bool_t fromFile)
{
   TClonesArray &list = *fList;
   TTVRecord *newrec = new(list[fRecords++])TTVRecord();
   if (!fromFile) newrec->FormFrom(fViewer);
   fCurrent = fRecords - 1;
   if (fRecords > 1) fViewer->ActivateButtons(kTRUE, kTRUE, kFALSE, kTRUE);
   else              fViewer->ActivateButtons(kTRUE, kFALSE, kFALSE, kTRUE);
   if (!fromFile) {
      TString name = "";
      if (strlen(newrec->GetZ())) name += newrec->GetZ();
      if (strlen(newrec->GetY())) {
         if (name.Length()) name += ":";
         name += newrec->GetY();
      }
      if (strlen(newrec->GetX())) {
         if (name.Length()) name += ":";
         name += newrec->GetX();
      }
      SetRecordName(name.Data());
   }
   return newrec;
}

////////////////////////////////////////////////////////////////////////////////
/// Return record at index i

TTVRecord *TTVSession::GetRecord(Int_t i)
{
   if (!fRecords) return 0;
   fCurrent = i;
   if (i < 0)           fCurrent = 0;
   if (i > fRecords-1)  fCurrent = fRecords - 1;
   if (fCurrent>0 && fCurrent<fRecords-1)
      fViewer->ActivateButtons(kTRUE, kTRUE, kTRUE, kTRUE);
   if (fCurrent == 0) {
      if (fRecords > 1) fViewer->ActivateButtons(kTRUE, kFALSE, kTRUE, kTRUE);
      else              fViewer->ActivateButtons(kTRUE, kFALSE, kFALSE, kTRUE);
   }
   if (fCurrent == fRecords-1) {
      if (fRecords > 1) fViewer->ActivateButtons(kTRUE, kTRUE, kFALSE, kTRUE);
      else              fViewer->ActivateButtons(kTRUE, kFALSE, kFALSE, kTRUE);
   }
   fViewer->SetCurrentRecord(fCurrent);
   return (TTVRecord *)fList->UncheckedAt(fCurrent);
}

////////////////////////////////////////////////////////////////////////////////
/// Set record name

void TTVSession::SetRecordName(const char *name)
{
   Int_t crt = fCurrent;
   TTVRecord *current = GetRecord(fCurrent);
   current->SetName(name);
   fViewer->UpdateCombo();
   fCurrent = crt;
   fViewer->SetCurrentRecord(fCurrent);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove current record from list

void TTVSession::RemoveLastRecord()
{
   if (!fRecords) return;
   TTVRecord *rec = (TTVRecord *)fList->UncheckedAt(fRecords);
   delete rec;
   fList->RemoveAt(fRecords--);
   if (fCurrent > fRecords-1) fCurrent = fRecords - 1;
   Int_t crt = fCurrent;
   fViewer->UpdateCombo();
   fCurrent = crt;
   if (!fRecords) {
      fViewer->ActivateButtons(kFALSE, kFALSE, kFALSE, kFALSE);
      return;
   }
   GetRecord(fCurrent);
}

////////////////////////////////////////////////////////////////////////////////
/// Display record rec

void TTVSession::Show(TTVRecord *rec)
{
   rec->PlugIn(fViewer);
   fViewer->ExecuteDraw();
   if (rec->HasUserCode() && rec->MustExecuteCode()) rec->ExecuteUserCode();
   fViewer->SetHistogramTitle(rec->GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Save the TTVSession in a C++ macro file

void TTVSession::SaveSource(std::ofstream &out)
{
   out<<"//--- session object"<<std::endl;
   out<<"   TTVSession* tv_session = new TTVSession(treeview);"<<std::endl;
   out<<"   treeview->SetSession(tv_session);"<<std::endl;
   TTVRecord *record;
   for (Int_t i=0; i<fRecords; i++) {
      record = GetRecord(i);
      record->SaveSource(out);
   }
   out<<"//--- Connect first record"<<std::endl;
   out<<"   tv_session->First();"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Updates current record according to new X, Y, Z settings

void TTVSession::UpdateRecord(const char *name)
{
   TTVRecord *current = (TTVRecord *)fList->UncheckedAt(fCurrent);
   current->FormFrom(fViewer);
   SetRecordName(name);
}

