/// \file ROOT/TFitPanel.cxx
/// \ingroup WebGui ROOT7
/// \author Sergey Linev <S.Linev@gsi.de>
/// \date 2017-10-24
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/TFitPanel.hxx"

#include <ROOT/TWebWindowsManager.hxx>
#include <ROOT/TLogger.hxx>
#include "ROOT/TDirectory.hxx"

#include "TString.h"
#include "TROOT.h"
#include "TBufferJSON.h"

std::shared_ptr<ROOT::Experimental::TWebWindow> ROOT::Experimental::TFitPanel::GetWindow()
{
   if (!fWindow) {
      fWindow = TWebWindowsManager::Instance()->CreateWindow(false);

      fWindow->SetPanelName("FitPanel");

      fWindow->SetDataCallBack(std::bind(&TFitPanel::ProcessData, this, std::placeholders::_1, std::placeholders::_2));
   }

   return fWindow;
}


void ROOT::Experimental::TFitPanel::Show(const std::string &where)
{
   GetWindow()->Show(where);
}

void ROOT::Experimental::TFitPanel::Hide()
{
   if (!fWindow) return;

   fWindow->CloseConnections();
}

void ROOT::Experimental::TFitPanel::ProcessData(unsigned connid, const std::string &arg)
{
   if (arg == "CONN_READY") {
      fConnId = connid;
      printf("Connection established %u\n", fConnId);
      fWindow->Send("INITDONE", fConnId);

      TFitPanelModel model;
      model.fDataNames.push_back(ComboBoxItem("1","RootData1"));
      model.fDataNames.push_back(ComboBoxItem("2","RootData2"));
      model.fDataNames.push_back(ComboBoxItem("3","RootData3"));
      model.fSelectDataId = "1";

      model.fModelNames.push_back(ComboBoxItem("1","RootModel1"));
      model.fModelNames.push_back(ComboBoxItem("2","RootModel2"));
      model.fModelNames.push_back(ComboBoxItem("3","RootModel3"));
      model.fSelectModelId = "3";

      TString json = TBufferJSON::ConvertToJSON(&model, gROOT->GetClass("ROOT::Experimental::TFitPanelModel"));

      fWindow->Send(std::string("MODEL:") + json.Data(), fConnId);

      return;
   }

   if (arg == "CONN_CLOSED") {
      printf("Connection closed\n");
      fConnId = 0;
      return;
   }

   if (arg.find("DOFIT:")==0) {
      TString exec;
      exec.Form("((ROOT::Experimental::TFitPanel *) %p)->DoFit(%s);", this, arg.c_str()+6);
      printf("Execute %s\n", exec.Data());
      gROOT->ProcessLine(exec.Data());
      return;
   }
}

void ROOT::Experimental::TFitPanel::UseCanvas(std::shared_ptr<TCanvas> &canv)
{
   if (fCanvas) {
      R__ERROR_HERE("ShowIn") << "FitPanel already bound to the canvas - change is not yet supported";
      return;
   }

   fCanvas = canv;
}


/// method called from the UI
void ROOT::Experimental::TFitPanel::DoFit(const std::string &dname, const std::string &mname)
{
   printf("DoFit %s %s\n", dname.c_str(), mname.c_str());

   bool first_time = false;

   if (!fCanvas) {
      fCanvas = Experimental::TCanvas::Create("FitPanel Canvas");
      first_time = true;
   }

   if (!fFitHist) {

      // Create the histogram.
      auto xaxis = std::make_shared<ROOT::Experimental::TAxisConfig>(10, 0., 10.);

      fFitHist = std::make_shared<ROOT::Experimental::TH1D>(*xaxis.get());

      // Fill a few points.
      fFitHist->Fill(5);
      fFitHist->Fill(6);
      fFitHist->Fill(6);
      fFitHist->Fill(7);

      fCanvas->Draw(fFitHist).SetLineColor(Experimental::TColor::kBlue);

      // workaround to keep histogram in the lists
      ROOT::Experimental::TDirectory::Heap().Add("fitaxis", xaxis);

      if (first_time) {
         fCanvas->Show();
         //fCanvas->Update();
      } else {
         fCanvas->Modified();
         //fCanvas->Update();
      }
   }
}
