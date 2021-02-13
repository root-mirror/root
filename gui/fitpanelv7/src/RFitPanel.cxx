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

#include <ROOT/RFitPanel.hxx>

#include <ROOT/RMakeUnique.hxx>
#include <ROOT/RLogger.hxx>

#include "Fit/BinData.h"
#include "Fit/Fitter.h"

#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TMultiGraph.h"
#include "THStack.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TFitResult.h"
#include "TList.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBufferJSON.h"
#include "Math/Minimizer.h"
#include "HFitInterface.h"
#include "TColor.h"

#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std::string_literals;

/** \class ROOT::Experimental::RFitPanel
\ingroup webdisplay

web-based FitPanel prototype.
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor

ROOT::Experimental::RFitPanel::FitRes::FitRes(const std::string &_objid, std::unique_ptr<TF1> &_func, TFitResultPtr &_res)
   : objid(_objid), res(_res)
{
   std::swap(func, _func);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

ROOT::Experimental::RFitPanel::FitRes::~FitRes()
{
   // to avoid dependency from TF1
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

ROOT::Experimental::RFitPanel::RFitPanel(const std::string &title)
{
   model().fTitle = title;

   GetFunctionsFromSystem();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

ROOT::Experimental::RFitPanel::~RFitPanel()
{
   // to avoid dependency from TF1
}

////////////////////////////////////////////////////////////////////////////////
/// Returns RWebWindow instance, used to display FitPanel

std::shared_ptr<ROOT::Experimental::RWebWindow> ROOT::Experimental::RFitPanel::GetWindow()
{
   if (!fWindow) {
      fWindow = RWebWindow::Create();

      fWindow->SetPanelName("rootui5.fitpanel.view.FitPanel");

      fWindow->SetCallBacks(
            [this](unsigned connid)
            {
               fConnId = connid;
               fWindow->Send(fConnId, "INITDONE");
               if (!model().fInitialized)
                  SelectObject("$$$");
               SendModel();
            },
            [this](unsigned connid, const std::string &arg) { ProcessData(connid, arg); },
            [this](unsigned) { fConnId = 0; });

      fWindow->SetGeometry(400, 650); // configure predefined geometry
   }

   return fWindow;
}

////////////////////////////////////////////////////////////////////////////////
/// Update list of available data

void ROOT::Experimental::RFitPanel::UpdateDataSet()
{
   auto &m = model();

   m.fDataSet.clear();

   for (auto &obj : fObjects)
      m.fDataSet.emplace_back("Panel", "panel::"s + obj->GetName(), Form("%s::%s", obj->ClassName(), obj->GetName()));

   if (gDirectory) {
      TIter iter(gDirectory->GetList());
      TObject *obj = nullptr;
      while ((obj = iter()) != nullptr) {
         if (GetFitObjectType(obj) != RFitPanelModel::kObjectNotSupported)
            m.fDataSet.emplace_back("gDirectory", "gdir::"s + obj->GetName(), Form("%s::%s", obj->ClassName(), obj->GetName()));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Select object for fitting

void ROOT::Experimental::RFitPanel::SelectObject(const std::string &objid)
{
   UpdateDataSet();

   auto &m = model();

   std::string id = objid;
   if (id.compare("$$$") == 0) {
      if (m.fDataSet.size() > 0)
         id = m.fDataSet[0].id;
      else
         id.clear();
   }

   TObject *obj = GetSelectedObject(id);
   auto kind = GetFitObjectType(obj);

   m.SetObjectKind(kind);

   TH1 *hist = nullptr;
   switch (kind) {
      case RFitPanelModel::kObjectHisto:
         hist = (TH1*)obj;
         break;

      case RFitPanelModel::kObjectGraph:
         hist = ((TGraph*)obj)->GetHistogram();
         break;

      case RFitPanelModel::kObjectMultiGraph:
         hist = ((TMultiGraph*)obj)->GetHistogram();
         break;

      case RFitPanelModel::kObjectGraph2D:
         hist = ((TGraph2D*)obj)->GetHistogram("empty");
         break;

      case RFitPanelModel::kObjectHStack:
         hist = (TH1 *)((THStack *)obj)->GetHists()->First();
         break;

      //case RFitPanelModel::kObjectTree:
      //   m.fFitMethods = {{ kFP_MUBIN, "Unbinned Likelihood" }};
      //   m.fFitMethod = kFP_MUBIN;
      //   break;

      default:
         break;
   }

   if (!obj)
      m.fSelectedData = "";
   else
      m.fSelectedData = id;

   m.fInitialized = true;

   // update list of data

   m.UpdateRange(hist);

   UpdateFunctionsList();

   std::string selfunc = m.fSelectedFunc;

   if (!m.HasFunction(selfunc)) {
      if (m.fFuncList.size() > 0)
         selfunc = m.fFuncList[0].id;
      else
         selfunc.clear();
   }

   SelectFunction(selfunc);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns object based on it string id
/// Searches either in gDirectory or in internal panel list

TObject *ROOT::Experimental::RFitPanel::GetSelectedObject(const std::string &objid)
{
   if (objid.compare(0,6,"gdir::") == 0) {
      std::string name = objid.substr(6);
      if (gDirectory)
         return gDirectory->GetList()->FindObject(name.c_str());
   } else if (objid.compare(0,7,"panel::") == 0) {
      std::string name = objid.substr(7);
      for (auto &item : fObjects)
         if (name.compare(item->GetName()) == 0)
            return item;
   }

   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kind of object

ROOT::Experimental::RFitPanelModel::EFitObjectType ROOT::Experimental::RFitPanel::GetFitObjectType(TObject *obj)
{
   if (!obj)
      return RFitPanelModel::kObjectNone;

   if (obj->InheritsFrom(TH1::Class()))
      return RFitPanelModel::kObjectHisto;

   if (obj->InheritsFrom(TGraph::Class()))
      return RFitPanelModel::kObjectGraph;

   if (obj->InheritsFrom(TGraph2D::Class()))
      return RFitPanelModel::kObjectGraph2D;

   if (obj->InheritsFrom(THStack::Class()))
      return RFitPanelModel::kObjectGraph2D;

   if (obj->InheritsFrom(TMultiGraph::Class()))
      return RFitPanelModel::kObjectGraph2D;

   return RFitPanelModel::kObjectNotSupported;
}

////////////////////////////////////////////////////////////////////////////////
/// Update list of available functions

void ROOT::Experimental::RFitPanel::UpdateFunctionsList()
{
   auto &m = model();

   m.fFuncList.clear();

   if (m.fDim == 1) {
      m.fFuncList = { {"gaus"}, {"gausn"}, {"expo"}, {"landau"}, {"landaun"},
                    {"pol0"},{"pol1"},{"pol2"},{"pol3"},{"pol4"},{"pol5"},{"pol6"},{"pol7"},{"pol8"},{"pol9"},
                    {"cheb0"}, {"cheb1"}, {"cheb2"}, {"cheb3"}, {"cheb4"}, {"cheb5"}, {"cheb6"}, {"cheb7"}, {"cheb8"}, {"cheb9"} };
   } else if (m.fDim == 2) {
      m.fFuncList = { {"xygaus"}, {"bigaus"}, {"xyexpo"}, {"xylandau"}, {"xylandaun"} };
   }

   for (auto &func : fSystemFuncs) {
      m.fFuncList.emplace_back("System", "system::"s + func->GetName(), func->GetName());
   }

   for (auto &entry : fPrevRes) {
      if (entry.objid == m.fSelectedData)
         m.fFuncList.emplace_back("Previous", "previous::"s + entry.func->GetName(), entry.func->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Select fit function

void ROOT::Experimental::RFitPanel::SelectFunction(const std::string &funcid)
{
   model().SelectedFunc(funcid, FindFunction(funcid));

   model().UpdateAdvanced(FindFitResult(funcid));
}

////////////////////////////////////////////////////////////////////////////////
/// Assign histogram to use with fit panel - without ownership

void ROOT::Experimental::RFitPanel::AssignHistogram(TH1 *hist)
{
   fObjects.emplace_back(hist);
   SelectObject("panel::"s + hist->GetName());
   SendModel();
}

/// Assign histogram name to use with fit panel - it should be available in gDirectory

void ROOT::Experimental::RFitPanel::AssignHistogram(const std::string &hname)
{
   SelectObject("gdir::" + hname);
   SendModel();
}

////////////////////////////////////////////////////////////////////////////////
/// assign canvas to use for drawing results of fitting or showing fitpanel itself

void ROOT::Experimental::RFitPanel::AssignCanvas(std::shared_ptr<RCanvas> &canv)
{
   if (!fCanvas) {
      fCanvas = canv;
   } else {
      R__LOG_ERROR(FitPanelLog()) << "FitPanel already bound to the canvas - change is not yet supported";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// assign histogram for fitting

void ROOT::Experimental::RFitPanel::AssignHistogram(std::shared_ptr<RH1D> &hist)
{
   fFitHist = hist;
}

////////////////////////////////////////////////////////////////////////////////
/// Show FitPanel

void ROOT::Experimental::RFitPanel::Show(const std::string &where)
{
   GetWindow()->Show(where);
}

////////////////////////////////////////////////////////////////////////////////
/// Hide FitPanel

void ROOT::Experimental::RFitPanel::Hide()
{
   if (fWindow)
      fWindow->CloseConnections();
}

////////////////////////////////////////////////////////////////////////////////
/// Return reference on model object
/// Model created if was not exists before


ROOT::Experimental::RFitPanelModel &ROOT::Experimental::RFitPanel::model()
{
   if (!fModel) {
      fModel = std::make_unique<RFitPanelModel>();
      fModel->Initialize();
   }

   return *fModel.get();
}

////////////////////////////////////////////////////////////////////////////////
/// Send model object to the client

void ROOT::Experimental::RFitPanel::SendModel()
{
   if (fWindow && (fConnId > 0)) {
      TString json = TBufferJSON::ToJSON(&model());
      fWindow->Send(fConnId, "MODEL:"s + json.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process data from FitPanel
/// OpenUI5-based FitPanel sends commands or status changes

void ROOT::Experimental::RFitPanel::ProcessData(unsigned, const std::string &arg)
{
   if (arg == "RELOAD") {

      GetFunctionsFromSystem();

      UpdateDataSet();
      UpdateFunctionsList();

      SendModel();

   } else if (arg.compare(0, 7, "UPDATE:") == 0) {

      if (UpdateModel(arg.substr(7)) > 0)
         SendModel();

   } else if (arg.compare(0, 6, "DOFIT:") == 0) {

      if (UpdateModel(arg.substr(6)) >= 0)
         if (DoFit())
            SendModel();

   } else if (arg.compare(0, 7, "DODRAW:") == 0) {

      if (UpdateModel(arg.substr(7)) >= 0)
         if (DoDraw())
            SendModel();

   } else if (arg.compare(0, 8, "SETPARS:") == 0) {

      auto info = TBufferJSON::FromJSON<RFitPanelModel::RFuncParsList>(arg.substr(8));

      if (info) {
         TF1 *func = FindFunction(info->id);

         // copy all parameters back to the function
         if (func)
            info->SetParameters(func);
      }

   }
}

////////////////////////////////////////////////////////////////////////////////
/// Search for existing functions, ownership still belongs to FitPanel or global lists

TF1 *ROOT::Experimental::RFitPanel::FindFunction(const std::string &id)
{
   if (id.compare(0,8,"system::") == 0) {
      std::string name = id.substr(8);

      for (auto &item : fSystemFuncs)
         if (name.compare(item->GetName()) == 0)
            return item.get();
   }

   if (id.compare(0,10,"previous::") == 0) {
      std::string name = id.substr(10);

      for (auto &entry : fPrevRes)
         if (name.compare(entry.func->GetName()) == 0)
            return entry.func.get();
   }

   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates new instance to make fitting

TFitResult *ROOT::Experimental::RFitPanel::FindFitResult(const std::string &id)
{
   if (id.compare(0,10,"previous::") == 0) {
      std::string name = id.substr(10);

      for (auto &entry : fPrevRes)
         if (name.compare(entry.func->GetName()) == 0)
            return entry.res.Get();
   }

   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates new instance to make fitting

std::unique_ptr<TF1> ROOT::Experimental::RFitPanel::GetFitFunction(const std::string &funcname)
{
   std::unique_ptr<TF1> res;

   TF1 *func = FindFunction(funcname);

   if (func) {
      // Now we make a copy.
      res.reset((TF1*) func->IsA()->New());
      func->Copy(*res);
   } else if (funcname.compare(0,6,"dflt::") == 0) {

      std::string formula = funcname.substr(6);

      ROOT::Fit::DataRange drange = model().GetRanges();

      double xmin, xmax, ymin, ymax, zmin, zmax;
      drange.GetRange(xmin, xmax, ymin, ymax, zmin, zmax);

      if ( model().fDim == 1 || model().fDim == 0 ) {
         res.reset(new TF1(formula.c_str(), formula.c_str(), xmin, xmax));
      } else if ( model().fDim == 2 ) {
         res.reset(new TF2(formula.c_str(), formula.c_str(), xmin, xmax, ymin, ymax));
      } else if ( model().fDim == 3 ) {
         res.reset(new TF3(formula.c_str(), formula.c_str(), xmin, xmax, ymin, ymax, zmin, zmax));
      }
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Update fit model
/// returns -1 if JSON fails
/// return 0 if nothing large changed
/// return 1 if important selection are changed and client need to be updated

int ROOT::Experimental::RFitPanel::UpdateModel(const std::string &json)
{
   auto m = TBufferJSON::FromJSON<RFitPanelModel>(json);

   if (!m) {
      R__LOG_ERROR(FitPanelLog()) << "Fail to parse JSON for RFitPanelModel";
      return -1;
   }

   m->fInitialized = true;

   int res = 0; // nothing changed

   if (model().fSelectedData != m->fSelectedData) {
      res |= 1;
   }

   if (model().fSelectedFunc != m->fSelectedFunc) {
      res |= 2;
   }

   std::swap(fModel, m); // finally replace model

   if (res & 1)
      SelectObject(model().fSelectedData);

   if (res != 0)
      SelectFunction(model().fSelectedFunc);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Copies f into a new TF1 to be stored in the fitpanel with it's
/// own ownership. This is taken from Fit::StoreAndDrawFitFunction in
/// HFitImpl.cxx

TF1* ROOT::Experimental::RFitPanel::copyTF1(TF1* f)
{
   double xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;

   // no need to use kNotGlobal bit. TF1::Copy does not add in the list by default
   if ( dynamic_cast<TF3*>(f)) {
      TF3* fnew = (TF3*)f->IsA()->New();
      f->Copy(*fnew);
      f->GetRange(xmin,ymin,zmin,xmax,ymax,zmax);
      fnew->SetRange(xmin,ymin,zmin,xmax,ymax,zmax);
      fnew->SetParent( nullptr );
      fnew->AddToGlobalList(false);
      return fnew;
   } else if ( dynamic_cast<TF2*>(f) != 0 ) {
      TF2* fnew = (TF2*)f->IsA()->New();
      f->Copy(*fnew);
      f->GetRange(xmin,ymin,xmax,ymax);
      fnew->SetRange(xmin,ymin,xmax,ymax);
      fnew->Save(xmin,xmax,ymin,ymax,0,0);
      fnew->SetParent( nullptr );
      fnew->AddToGlobalList(false);
      return fnew;
   }

   TF1* fnew = (TF1*)f->IsA()->New();
   f->Copy(*fnew);
   f->GetRange(xmin,xmax);
   fnew->SetRange(xmin,xmax);
   // This next line is added, as fnew-Save fails with gausND! As
   // the number of dimensions is unknown...
   if ( '\0' != fnew->GetExpFormula()[0] )
      fnew->Save(xmin,xmax,0,0,0,0);
   fnew->SetParent( nullptr );
   fnew->AddToGlobalList(false);
   return fnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Looks for all the functions registered in the current ROOT
/// session.

void ROOT::Experimental::RFitPanel::GetFunctionsFromSystem()
{

   fSystemFuncs.clear();

   // Be carefull not to store functions that will be in the
   // predefined section
   std::vector<std::string> fnames = { "gaus" ,   "gausn", "expo", "landau",
                                       "landaun", "pol0",  "pol1", "pol2",
                                       "pol3",    "pol4",  "pol5", "pol6",
                                       "pol7",    "pol8",  "pol9", "user" };

   // No go through all the objects registered in gROOT
   TIter functionsIter(gROOT->GetListOfFunctions());
   TObject* obj;
   while( (obj = functionsIter()) != nullptr ) {
      // And if they are TF1s
      if ( TF1* func = dynamic_cast<TF1*>(obj) ) {
         bool addFunction = true;
         // And they are not already registered in fSystemFunc
         for ( auto &name : fnames) {
            if ( name.compare(func->GetName()) == 0 ) {
               addFunction = false;
               break;
            }
         }
         // Add them.
         if ( addFunction )
            fSystemFuncs.emplace_back( copyTF1(func) );
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns pad where histogram is drawn
/// If canvas not exists, create new one

TPad *ROOT::Experimental::RFitPanel::GetDrawPad(TObject *obj, bool force)
{
   if (!obj || (!force && (model().fNoDrawing || model().fNoStoreDraw)))
      return nullptr;

   std::function<TPad*(TPad*)> check = [&](TPad *pad) {
      TPad *res = nullptr;
      if (!pad) return res;
      if (!fPadName.empty() && (fPadName.compare(pad->GetName()) == 0)) return pad;
      TIter next(pad->GetListOfPrimitives());
      TObject *prim = nullptr;
      while (!res && (prim = next())) {
         res = (prim == obj) ? pad : check(dynamic_cast<TPad *>(prim));
      }
      return res;
   };

   if (!fCanvName.empty()) {
      auto drawcanv = dynamic_cast<TCanvas *> (gROOT->GetListOfCanvases()->FindObject(fCanvName.c_str()));
      auto drawpad = check(drawcanv);
      if (drawpad) {
         drawpad->cd();
         return drawpad;
      }
      if (drawcanv) {
         drawcanv->Clear();
         drawcanv->cd();
         obj->Draw();
         return drawcanv;
      }
      fCanvName.clear();
      fPadName.clear();
   }

   TObject *c = nullptr;
   TIter nextc(gROOT->GetListOfCanvases());
   while ((c = nextc())) {
      auto drawpad = check(dynamic_cast<TCanvas* >(c));
      if (drawpad) {
         drawpad->cd();
         fCanvName = c->GetName();
         fPadName = drawpad->GetName();
         return drawpad;
      }
   }

   auto canv = gROOT->MakeDefCanvas();
   canv->SetName("fpc");
   canv->SetTitle("Fit panel drawings");
   fPadName = fCanvName = canv->GetName();

   canv->cd();
   obj->Draw();

   return canv;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform fitting using current model settings
/// Returns true if any action was done

bool ROOT::Experimental::RFitPanel::DoFit()
{
   auto &m = model();

   TObject *obj = GetSelectedObject(m.fSelectedData);
   if (!obj) return false;
   auto kind = GetFitObjectType(obj);

   auto f1 = GetFitFunction(m.fSelectedFunc);
   if (!f1) return false;

   auto drange = m.GetRanges();
   auto minOption = m.GetMinimizerOptions();
   auto fitOpts = m.GetFitOptions();
   auto drawOpts = m.GetDrawOption();

   fitOpts.StoreResult = 1;

   TVirtualPad *save = gPad;

   auto pad = GetDrawPad(obj);

   TFitResultPtr res;

   switch (kind) {
      case RFitPanelModel::kObjectHisto: {

         TH1 *hist = dynamic_cast<TH1*>(obj);
         if (hist)
            res = ROOT::Fit::FitObject(hist, f1.get(), fitOpts, minOption, drawOpts, drange);

         break;
      }
      case RFitPanelModel::kObjectGraph: {

         TGraph *gr = dynamic_cast<TGraph*>(obj);
         if (gr)
            res = ROOT::Fit::FitObject(gr, f1.get(), fitOpts, minOption, drawOpts, drange);
         break;
      }
      case RFitPanelModel::kObjectMultiGraph: {

         TMultiGraph *mg = dynamic_cast<TMultiGraph*>(obj);
         if (mg)
            res = ROOT::Fit::FitObject(mg, f1.get(), fitOpts, minOption, drawOpts, drange);

         break;
      }
      case RFitPanelModel::kObjectGraph2D: {

         TGraph2D *g2d = dynamic_cast<TGraph2D*>(obj);
         if (g2d)
            res = ROOT::Fit::FitObject(g2d, f1.get(), fitOpts, minOption, drawOpts, drange);

         break;
      }
      case RFitPanelModel::kObjectHStack: {
         // N/A
         break;
      }

      default: {
         // N/A
         break;
      }
   }

   // After fitting function appears in global list
   if (f1 && gROOT->GetListOfFunctions()->FindObject(f1.get()))
      gROOT->GetListOfFunctions()->Remove(f1.get());

   if (m.fSame && f1 && pad) {
      TF1 *copy = copyTF1(f1.get());
      copy->SetBit(kCanDelete);
      copy->Draw("same");
   }

   if (pad) {
      pad->Modified();
      pad->Update();
   }

   std::string funcname = f1->GetName();
   if ((funcname.compare(0,4,"prev") == 0) && (funcname.find("-") > 4))
      funcname.erase(0, funcname.find("-") + 1);
   funcname = "prev"s + std::to_string(fPrevRes.size() + 1) + "-"s + funcname;
   f1->SetName(funcname.c_str());

   fPrevRes.emplace_back(m.fSelectedData, f1, res);

   UpdateFunctionsList();

   SelectFunction("previous::"s + funcname);

   if (save && (gPad != save))
      gPad = save;

   return true; // provide client with latest changes
}

////////////////////////////////////////////////////////////////////////////////
/// Extract color from string
/// Should be coded as #ff00ff string

Color_t ROOT::Experimental::RFitPanel::GetColor(const std::string &colorid)
{
   if ((colorid.length() != 7) || (colorid.compare(0,1,"#") != 0)) return 0;

   return TColor::GetColor(colorid.c_str());
}

////////////////////////////////////////////////////////////////////////////////
/// Create confidence levels drawing
/// tab. Then it call Virtual Fitter to perform it.

TObject *ROOT::Experimental::RFitPanel::MakeConfidenceLevels(TFitResult *result)
{
   if (!result)
      return nullptr;

   // try to use provided method
   // auto conf = result->GetConfidenceIntervals();
   // printf("GET INTERVALS %d\n", (int) conf.size());

   const auto *function = result->FittedFunction();
   if (!function) {
      R__LOG_ERROR(FitPanelLog()) << "Fit Function does not exist!";
      return nullptr;
   }

   const auto *data = result->FittedBinData();
   if (!data) {
      R__LOG_ERROR(FitPanelLog()) << "Unbinned data set cannot draw confidence levels.";
      return nullptr;
   }

   std::vector<Double_t> ci(data->Size());
   result->GetConfidenceIntervals(*data, &ci[0], model().fConfidenceLevel);

   if (model().fDim == 1) {
      TGraphErrors *g = new TGraphErrors(ci.size());
      for (unsigned int i = 0; i < ci.size(); ++i) {
         const Double_t *x = data->Coords(i);
         const Double_t y = (*function)(x);
         g->SetPoint(i, *x, y);
         g->SetPointError(i, 0, ci[i]);
      }
      // std::ostringstream os;
      // os << "Confidence Intervals with " << conflvl << " conf. band.";
      // g->SetTitle(os.str().c_str());
      g->SetTitle("Confidence Intervals with");

      auto icol = GetColor(model().fConfidenceColor);
      g->SetLineColor(icol);
      g->SetFillColor(icol);
      g->SetFillStyle(3001);
      return g;
   } else if (model().fDim == 2) {
      TGraph2DErrors *g = new TGraph2DErrors(ci.size());
      for (unsigned int i = 0; i < ci.size(); ++i) {
         const Double_t *x = data->Coords(i);
         const Double_t y = (*function)(x);
         g->SetPoint(i, x[0], x[1], y);
         g->SetPointError(i, 0, 0, ci[i]);
      }
      // std::ostringstream os;
      // os << "Confidence Intervals with " << fConfLevel->GetNumber() << " conf. band.";
      // g->SetTitle(os.str().c_str());

      g->SetTitle("Confidence Intervals with");

      auto icol = GetColor(model().fConfidenceColor);
      g->SetLineColor(icol);
      g->SetFillColor(icol);
      g->SetFillStyle(3001);
      return g;
   }

   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform drawing using current model settings
/// Returns true if any action was done

bool ROOT::Experimental::RFitPanel::DoDraw()
{
   auto &m = model();

   TObject *obj = GetSelectedObject(m.fSelectedData);
   if (!obj) return false;

   TObject *drawobj = nullptr;
   std::string drawopt;
   bool superimpose = true, objowner = true;

   if (m.fHasAdvanced && (m.fSelectedTab == "Advanced")) {

      TFitResult *res = FindFitResult(m.fSelectedFunc);
      if (!res) return false;

      if (m.fAdvancedTab == "Contour") {

         superimpose = m.fContourSuperImpose;
         int par1 = std::stoi(m.fContourPar1Id);
         int par2 = std::stoi(m.fContourPar2Id);

         TGraph *graph = new TGraph(m.fContourPoints);

         if (!res->Contour(par1, par2, graph, m.fConfidenceLevel)) {
            delete graph;
            return false;
         }

         auto fillcolor = GetColor(m.fContourColor);
         graph->SetFillColor(fillcolor);
         graph->GetXaxis()->SetTitle( res->ParName(par1).c_str() );
         graph->GetYaxis()->SetTitle( res->ParName(par2).c_str() );

         drawobj = graph;
         drawopt = superimpose ? "LF" : "ALF";

      } else if (m.fAdvancedTab == "Scan") {

         int par = std::stoi(m.fScanId);
         TGraph *graph = new TGraph( m.fScanPoints);
         if (!res->Scan( par, graph, m.fScanMin, m.fScanMax)) {
            delete graph;
            return false;
         }
         auto linecolor = GetColor(m.fScanColor);
         if (!linecolor) linecolor = kBlue;
         graph->SetLineColor(linecolor);
         graph->SetLineWidth(2);
         graph->GetXaxis()->SetTitle(res->ParName(par).c_str());
         graph->GetYaxis()->SetTitle("FCN" );

         superimpose = false;
         drawobj = graph;
         drawopt = "ALF";

      } else if (m.fAdvancedTab == "Confidence") {

         drawobj = MakeConfidenceLevels(res);
         drawopt = "C3same";

      } else {
         return false;
      }
   } else {

       // find already existing functions, not try to create something new
       TF1 *func = FindFunction(m.fSelectedFunc);

       // when "Pars" tab is selected, automatically update function parameters
       if (func && (m.fSelectedTab.compare("Pars") == 0) && (m.fSelectedFunc == m.fFuncPars.id))
           m.fFuncPars.SetParameters(func);

       drawobj = func;
       drawopt = "same";
       objowner = true;
   }

   if (!drawobj)
      return false;

   auto pad = GetDrawPad(obj, true);
   if (!pad) {
      if (objowner)
         delete drawobj;
      return false;
   }

   if (!superimpose)
      pad->Clear();

   if (objowner)
      drawobj->SetBit(kCanDelete);

   drawobj->Draw(drawopt.c_str());

   pad->Modified();
   pad->Update();

   return true;
}
