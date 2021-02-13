// Authors: Sergey Linev <S.Linev@gsi.de> Iliana Betsou <Iliana.Betsou@cern.ch>
// Date: 2019-04-11
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RFitPanel
#define ROOT_RFitPanel

#include <ROOT/RWebWindow.hxx>

#include <ROOT/RFitPanelModel.hxx>

#include <ROOT/RCanvas.hxx>

#include "ROOT/RHist.hxx"

#include "TFitResultPtr.h"

#include <memory>
#include <vector>
#include <list>

class TPad;
class TH1;
class TF1;

namespace ROOT {
namespace Experimental {

class RFitPanel {

   std::unique_ptr<RFitPanelModel> fModel;

   std::vector<TObject*> fObjects;    ///<! objects provided directly to panel for fitting
   std::string fCanvName;             ///<! v6 canvas name used to display fit, will be created if not exists
   std::string fPadName;              ///<! v6 pad name in the canvas, where object is (was) drawn

   std::shared_ptr<RCanvas> fCanvas; ///<! v7 canvas used to display results
   std::shared_ptr<RH1D> fFitHist;   ///<! v7 histogram for fitting

   std::shared_ptr<RWebWindow> fWindow;  ///<! configured display
   unsigned fConnId{0};                  ///<! client connection id

   std::vector<std::unique_ptr<TF1>> fSystemFuncs; ///<! local copy of all internal system funcs

   struct FitRes {
      std::string objid;         // object used for fitting
      std::unique_ptr<TF1> func; // function
      TFitResultPtr res;          // result
      FitRes() = default;
      FitRes(const std::string &_objid, std::unique_ptr<TF1> &_func, TFitResultPtr &_res);
      ~FitRes();
   };

   std::list<FitRes> fPrevRes; ///<! all previous functions used for fitting

   TF1 *copyTF1(TF1 *f);

   void GetFunctionsFromSystem();

   void ProcessData(unsigned connid, const std::string &arg);

   RFitPanelModel &model();

   int UpdateModel(const std::string &json);

   void SelectObject(const std::string &objid);
   TObject *GetSelectedObject(const std::string &objid);
   RFitPanelModel::EFitObjectType GetFitObjectType(TObject *obj);
   void UpdateDataSet();

   void UpdateFunctionsList();
   TF1 *FindFunction(const std::string &funcid);
   TFitResult *FindFitResult(const std::string &funcid);
   std::unique_ptr<TF1> GetFitFunction(const std::string &funcid);
   void SelectFunction(const std::string &funcid);

   TObject *MakeConfidenceLevels(TFitResult *res);

   Color_t GetColor(const std::string &colorid);

   TPad *GetDrawPad(TObject *obj, bool force = false);

   void SendModel();

   bool DoFit();
   bool DoDraw();

public:
   RFitPanel(const std::string &title = "Fit panel");
   ~RFitPanel();

   // method required when any panel want to be inserted into the RCanvas
   std::shared_ptr<RWebWindow> GetWindow();

   void AssignHistogram(TH1 *hist);

   void AssignHistogram(const std::string &hname);

   void AssignCanvas(const std::string &cname) { fCanvName = cname; }

   void AssignCanvas(std::shared_ptr<RCanvas> &canv);

   void AssignHistogram(std::shared_ptr<RH1D> &hist);

   /// show FitPanel in specified place
   void Show(const std::string &where = "");

   /// hide FitPanel
   void Hide();
};

} // namespace Experimental
} // namespace ROOT

#endif
