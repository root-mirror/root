// Author: Sergey Linev, GSI  14/01/2021

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef RGeomViewerWidget_H
#define RGeomViewerWidget_H

#include <QWidget>
#include <QWebEngineView>

#include <memory>

namespace ROOT {
namespace Experimental {
class REveGeomViewer;
}
}

class RGeomViewerWidget : public QWidget {

   Q_OBJECT

protected:

   void resizeEvent(QResizeEvent *event) override;

   QWebEngineView *fView{nullptr};  ///< qt webwidget to show

   std::shared_ptr<ROOT::Experimental::REveGeomViewer> fGeomViewer;

public:
   RGeomViewerWidget(QWidget *parent = nullptr);
   virtual ~RGeomViewerWidget();

   /// returns geometry viewer
   auto getGeomViewer() { return fGeomViewer; }
};

#endif
