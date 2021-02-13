// Author: Sergey Linev, GSI  13/01/2021

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TCanvasWidget.h"

#include "TCanvas.h"
#include "TROOT.h"
#include "TClass.h"
#include "TEnv.h"

#include "TWebCanvas.h"

TCanvasWidget::TCanvasWidget(QWidget *parent) : QWidget(parent)
{
   setObjectName( "TCanvasWidget");

   setSizeIncrement( QSize( 100, 100 ) );

   setUpdatesEnabled( true );
   setMouseTracking(true);

   setFocusPolicy( Qt::TabFocus );
   setCursor( Qt::CrossCursor );

   setAcceptDrops(true);

   static int wincnt = 1;

   fCanvas = new TCanvas(kFALSE);
   fCanvas->SetName(Form("Canvas%d", wincnt++));
   fCanvas->SetTitle("title");
   fCanvas->ResetBit(TCanvas::kShowEditor);
   fCanvas->SetCanvas(fCanvas);
   fCanvas->SetBatch(kTRUE); // mark canvas as batch

   gPad = fCanvas;

   bool read_only = (gEnv->GetValue("WebGui.FullCanvas", (Int_t) 0) == 0);

   TWebCanvas *web = new TWebCanvas(fCanvas, "title", 0, 0, 800, 600, read_only);

   fCanvas->SetCanvasImp(web);

   SetPrivateCanvasFields(true);

   web->SetCanCreateObjects(kFALSE); // not yet create objects on server side

   web->SetUpdatedHandler([=]() { emit CanvasUpdated(); });

   web->SetActivePadChangedHandler([=](TPad *pad){ emit SelectedPadChanged(pad); });

   web->SetPadClickedHandler([=](TPad *pad, int x, int y) { emit PadClicked(pad,x,y); });

   web->SetPadDblClickedHandler([this](TPad *pad, int x, int y) { emit PadDblClicked(pad,x,y); });

   auto where = ROOT::Experimental::RWebDisplayArgs::GetQt5EmbedQualifier(this, "noopenui");

   web->ShowWebWindow(where);

   fView = findChild<QWebEngineView*>("RootWebView");
   if (!fView) {
      printf("FAIL TO FIND QWebEngineView - ROOT Qt5Web plugin does not work properly !!!!!\n");
      exit(11);
   }

   fView->resize(width(), height());
   fCanvas->SetCanvasSize(width(), height());
}

TCanvasWidget::~TCanvasWidget()
{
   if (fCanvas) {
      SetPrivateCanvasFields(false);

      gROOT->GetListOfCanvases()->Remove(fCanvas);

      fCanvas->Close();
      delete fCanvas;
      fCanvas = nullptr;
   }
}

void TCanvasWidget::SetPrivateCanvasFields(bool on_init)
{
   Long_t offset = TCanvas::Class()->GetDataMemberOffset("fCanvasID");
   if (offset > 0) {
      Int_t *id = (Int_t *)((char*) fCanvas + offset);
      if (*id == fCanvas->GetCanvasID()) *id = on_init ? 111222333 : -1;
   } else {
      printf("ERROR: Cannot modify fCanvasID data member\n");
   }

   offset = TCanvas::Class()->GetDataMemberOffset("fMother");
   if (offset > 0) {
      TPad **moth = (TPad **)((char*) fCanvas + offset);
      if (*moth == fCanvas->GetMother()) *moth = on_init ? fCanvas : nullptr;
   } else {
      printf("ERROR: Cannot set fMother data member in canvas\n");
   }
}

void TCanvasWidget::resizeEvent(QResizeEvent *event)
{
   fView->resize(width(), height());
   fCanvas->SetCanvasSize(width(), height());
}

void TCanvasWidget::activateEditor(TPad *pad, TObject *obj)
{
   TWebCanvas *cimp = dynamic_cast<TWebCanvas*> (fCanvas->GetCanvasImp());
   if (cimp) {
      cimp->ShowEditor(kTRUE);
      cimp->ActivateInEditor(pad, obj);
   }
}

void TCanvasWidget::setEditorVisible(bool flag)
{
   TCanvasImp *cimp = fCanvas->GetCanvasImp();
   if (cimp) cimp->ShowEditor(flag);
}

void TCanvasWidget::activateStatusLine()
{
   TCanvasImp *cimp = fCanvas->GetCanvasImp();
   if (cimp) cimp->ShowStatusBar(kTRUE);
}
