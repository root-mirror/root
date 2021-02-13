// Author: Sergey Linev, GSI  13/01/2021

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "RCanvasWidget.h"

#include <ROOT/RCanvas.hxx>
#include <ROOT/RWebDisplayArgs.hxx>

RCanvasWidget::RCanvasWidget(QWidget *parent) : QWidget(parent)
{

   setObjectName( "RCanvasWidget");

   setSizeIncrement( QSize( 100, 100 ) );

   setUpdatesEnabled( true );
   setMouseTracking(true);

   setFocusPolicy( Qt::TabFocus );
   setCursor( Qt::CrossCursor );

   setAcceptDrops(true);

   fCanvas = ROOT::Experimental::RCanvas::Create("ExampleCanvas");

   auto where = ROOT::Experimental::RWebDisplayArgs::GetQt5EmbedQualifier(this, "noopenui");

   fCanvas->Show(where);

   fView = findChild<QWebEngineView*>("RootWebView");
   if (!fView) {
      printf("FAIL TO FIND QWebEngineView - ROOT Qt5Web plugin does not work properly !!!!!\n");
      exit(11);
   }

   fView->resize(width(), height());
   fCanvas->SetSize({width(), height()});
}

RCanvasWidget::~RCanvasWidget()
{
   // remove canvas from global lists
   fCanvas->Remove();
}

void RCanvasWidget::resizeEvent(QResizeEvent *event)
{
   fView->resize(width(), height());
   fCanvas->SetSize({width(), height()});
}
