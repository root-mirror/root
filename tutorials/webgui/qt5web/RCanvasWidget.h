// Author: Sergey Linev, GSI  13/01/2021

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef RCanvasWidget_H
#define RCanvasWidget_H

#include <QWidget>
#include <QWebEngineView>

#include <memory>

namespace ROOT {
namespace Experimental {
class RCanvas;
}
}

class RCanvasWidget : public QWidget {

   Q_OBJECT

protected:

   void resizeEvent(QResizeEvent *event) override;

   QWebEngineView *fView{nullptr};  ///< qt webwidget to show

   std::shared_ptr<ROOT::Experimental::RCanvas> fCanvas;

public:
   RCanvasWidget(QWidget *parent = nullptr);
   virtual ~RCanvasWidget();

   /// returns canvas shown in the widget
   auto getCanvas() { return fCanvas; }


};

#endif
