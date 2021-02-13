// Author: Sergey Linev, GSI  13/01/2021

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveGeomViewer.hxx>
#include <ROOT/RWebDisplayArgs.hxx>
#include "RGeomViewerWidget.h"

RGeomViewerWidget::RGeomViewerWidget(QWidget *parent) : QWidget(parent)
{
   setObjectName( "RGeomViewerWidget");

   setSizeIncrement( QSize( 100, 100 ) );

   setUpdatesEnabled( true );
   setMouseTracking(true);

   setFocusPolicy( Qt::TabFocus );
   setCursor( Qt::CrossCursor );

   setAcceptDrops(true);

   fGeomViewer = std::make_shared<ROOT::Experimental::REveGeomViewer>();

   fGeomViewer->SetShowHierarchy(false);

   auto where = ROOT::Experimental::RWebDisplayArgs::GetQt5EmbedQualifier(this);

   fGeomViewer->Show(where);

   fView = findChild<QWebEngineView*>("RootWebView");
   if (!fView) {
      printf("FAIL TO FIND QWebEngineView - ROOT Qt5Web plugin does not work properly !!!!!\n");
      exit(11);
   }

   fView->resize(width(), height());
}

RGeomViewerWidget::~RGeomViewerWidget()
{
}

void RGeomViewerWidget::resizeEvent(QResizeEvent *event)
{
   fView->resize(width(), height());
}
