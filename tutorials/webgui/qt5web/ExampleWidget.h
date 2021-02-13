// Author: Sergey Linev, GSI  13/01/2021

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ExampleWidget_h
#define ExampleWidget_h

#include <QWidget>
#include "ui_ExampleWidget.h"

#include <memory>

class TH1F;
class TH2I;

class ExampleWidget : public QWidget, public Ui::ExampleWidget
{
   Q_OBJECT

   protected:

      TH1F *fHisto{nullptr};  ///< histogram for display in TCanvas
      std::shared_ptr<TH2I> fHisto2; ///< histogram for display in RCanvas

      void ImportCmsGeometry();

      void CreateDummyGeometry();

   public:

      ExampleWidget(QWidget *parent = nullptr, const char* name = nullptr);

      virtual ~ExampleWidget();

   public slots:

      void InfoButton_clicked();
      void CmsButton_clicked();
      void ExitButton_clicked();

};

#endif
