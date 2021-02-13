// Authors: Bertrand Bellenot <bertrand.bellenot@cern.ch> Sergey Linev <S.Linev@gsi.de>
// Date: 2019-02-28
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RBrowser.hxx>

#include <ROOT/Browsable/RGroup.hxx>
#include <ROOT/Browsable/RWrapper.hxx>
#include <ROOT/Browsable/RProvider.hxx>
#include <ROOT/Browsable/TObjectHolder.hxx>
#include <ROOT/Browsable/RSysFile.hxx>

#include <ROOT/RLogger.hxx>
#include <ROOT/RMakeUnique.hxx>
#include <ROOT/RFileDialog.hxx>

#include "RBrowserWidget.hxx"

#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TBufferJSON.h"
#include "TApplication.h"
#include "TRint.h"
#include "Getline.h"

#include <sstream>
#include <iostream>
#include <algorithm>
#include <memory>
#include <mutex>
#include <thread>
#include <fstream>

using namespace std::string_literals;

using namespace ROOT::Experimental;


class RBrowserEditorWidget : public RBrowserWidget {
public:

   bool fIsEditor{true};   ///<! either editor or image viewer
   std::string fTitle;
   std::string fFileName;
   std::string fContent;
   bool fFirstSend{false};  ///<! if editor content was send at least one
   std::string fItemPath;   ///<! item path in the browser

   RBrowserEditorWidget(const std::string &name, bool is_editor = true) : RBrowserWidget(name), fIsEditor(is_editor) {}
   virtual ~RBrowserEditorWidget() = default;

   void ResetConn() override { fFirstSend = false; }

   std::string GetKind() const override { return fIsEditor ? "editor"s : "image"s; }
   std::string GetTitle() override { return fTitle; }
   std::string GetUrl() override { return ""s; }

   void Show(const std::string &) override {}

   bool DrawElement(std::shared_ptr<Browsable::RElement> &elem, const std::string &) override
   {
      if (fIsEditor && elem->IsCapable(Browsable::RElement::kActEdit)) {
         auto code = elem->GetContent("text");
         if (!code.empty()) {
            fFirstSend = false;
            fContent = code;
            fTitle = elem->GetName();
            fFileName = elem->GetContent("filename");
         } else {
            auto json = elem->GetContent("json");
            if (!json.empty()) {
               fFirstSend = false;
               fContent = json;
               fTitle = elem->GetName() + ".json";
               fFileName = "";
            }
         }
         if (!fContent.empty()) {
            // page->fItemPath = item_path;
            return true;
         }
      }

      if (!fIsEditor && elem->IsCapable(Browsable::RElement::kActImage)) {
         auto img = elem->GetContent("image64");
         if (!img.empty()) {
            fFirstSend = false;
            fContent = img;
            fTitle = elem->GetName();
            fFileName = elem->GetContent("filename");
            // fItemPath = item_path;

            return true;
         }
      }

      return false;
   }

   std::string SendWidgetContent() override
   {
      if (fFirstSend) return ""s;

      fFirstSend = true;
      std::vector<std::string> args = { GetName(), fTitle, fFileName, fContent };

      std::string msg = fIsEditor ? "EDITOR:"s : "IMAGE:"s;
      msg += TBufferJSON::ToJSON(&args).Data();
      return msg;
   }

};


/** \class ROOT::Experimental::RBrowser
\ingroup rbrowser
\brief Web-based %ROOT file browser
*/

//////////////////////////////////////////////////////////////////////////////////////////////
/// constructor

RBrowser::RBrowser(bool use_rcanvas)
{
   SetUseRCanvas(use_rcanvas);

   auto comp = std::make_shared<Browsable::RGroup>("top","Root browser");

   auto seldir = Browsable::RSysFile::ProvideTopEntries(comp);

   std::unique_ptr<Browsable::RHolder> rootfold = std::make_unique<Browsable::TObjectHolder>(gROOT->GetRootFolder(), kFALSE);
   auto elem_root = Browsable::RProvider::Browse(rootfold);
   if (elem_root)
      comp->Add(std::make_shared<Browsable::RWrapper>("root", elem_root));

   std::unique_ptr<Browsable::RHolder> rootfiles = std::make_unique<Browsable::TObjectHolder>(gROOT->GetListOfFiles(), kFALSE);
   auto elem_files = Browsable::RProvider::Browse(rootfiles);
   if (elem_files)
      comp->Add(std::make_shared<Browsable::RWrapper>("ROOT Files", elem_files));

   fBrowsable.SetTopElement(comp);

   fBrowsable.SetWorkingPath(seldir);

   fWebWindow = RWebWindow::Create();
   fWebWindow->SetDefaultPage("file:rootui5sys/browser/browser.html");

   // this is call-back, invoked when message received via websocket
   fWebWindow->SetCallBacks([this](unsigned connid) { fConnId = connid; SendInitMsg(connid); },
                            [this](unsigned connid, const std::string &arg) { ProcessMsg(connid, arg); });
   fWebWindow->SetGeometry(1200, 700); // configure predefined window geometry
   fWebWindow->SetConnLimit(1); // the only connection is allowed
   fWebWindow->SetMaxQueueLength(30); // number of allowed entries in the window queue

   Show();

   // add first canvas by default

   //if (GetUseRCanvas())
   //   AddWidget("rcanvas");
   //else
   //   AddWidget("tcanvas");

   // AddWidget("geom");  // add geometry viewer at the beginning

   // AddWidget("editor"); // one can add empty editor if necessary
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// destructor

RBrowser::~RBrowser()
{
}


//////////////////////////////////////////////////////////////////////////////////////////////
/// Process browser request

std::string RBrowser::ProcessBrowserRequest(const std::string &msg)
{
   std::unique_ptr<RBrowserRequest> request;

   if (msg.empty()) {
      request = std::make_unique<RBrowserRequest>();
      request->first = 0;
      request->number = 100;
   } else {
      request = TBufferJSON::FromJSON<RBrowserRequest>(msg);
   }

   if (!request)
      return ""s;

   return "BREPL:"s + fBrowsable.ProcessRequest(*request.get());
}

/////////////////////////////////////////////////////////////////////////////////
/// Process file save command in the editor

void RBrowser::ProcessSaveFile(const std::string &fname, const std::string &content)
{
   if (fname.empty()) return;
   R__LOG_DEBUG(0, BrowserLog()) << "SaveFile " << fname << "  content length " << content.length();
   std::ofstream f(fname);
   f << content;
}

/////////////////////////////////////////////////////////////////////////////////
/// Process file save command in the editor

long RBrowser::ProcessRunMacro(const std::string &file_path)
{
   return gInterpreter->ExecuteMacro(file_path.c_str());
}

/////////////////////////////////////////////////////////////////////////////////
/// Process dbl click on browser item

std::string RBrowser::ProcessDblClick(std::vector<std::string> &args)
{
   args.pop_back(); // remove exec string, not used now

   std::string drawingOptions = args.back();
   args.pop_back(); // remove draw option

   auto path = fBrowsable.GetWorkingPath();
   path.insert(path.end(), args.begin(), args.end());

   R__LOG_DEBUG(0, BrowserLog()) << "DoubleClick " << Browsable::RElement::GetPathAsString(path);

   auto elem = fBrowsable.GetSubElement(path);
   if (!elem) return ""s;

   auto widget = GetActiveWidget();
   if (widget && widget->DrawElement(elem, drawingOptions)) {
      widget->SetPath(path);
      return widget->SendWidgetContent();
   }

   // check if element was drawn in other widget and just activate that widget
   auto iter = std::find_if(fWidgets.begin(), fWidgets.end(),
         [path](const std::shared_ptr<RBrowserWidget> &wg) { return path == wg->GetPath(); });

   if (iter != fWidgets.end())
      return "SELECT_WIDGET:"s + (*iter)->GetName();

   auto dflt_action = elem->GetDefaultAction();

   // check if object can be drawn in RCanvas even when default action is drawing in TCanvas
   if ((dflt_action == Browsable::RElement::kActDraw6) && GetUseRCanvas() && elem->IsCapable(Browsable::RElement::kActDraw7))
      dflt_action = Browsable::RElement::kActDraw7;

   std::string widget_kind;
   switch(dflt_action) {
      case Browsable::RElement::kActGeom: widget_kind = "geom"; break;
      case Browsable::RElement::kActDraw6: widget_kind = "tcanvas"; break;
      case Browsable::RElement::kActDraw7: widget_kind = "rcanvas"; break;
      case Browsable::RElement::kActEdit: widget_kind = "editor"; break;
      case Browsable::RElement::kActImage: widget_kind = "image"; break;
      default: widget_kind.clear();
   }

   if (!widget_kind.empty()) {
      auto new_widget = AddWidget(widget_kind);
      if (new_widget) {
         // draw object before client side is created - should not be a problem
         // after widget add in browser, connection will be established and data provided
         if (new_widget->DrawElement(elem, drawingOptions))
            new_widget->SetPath(path);
         return NewWidgetMsg(new_widget);
      }
   }

   if (elem->IsCapable(Browsable::RElement::kActBrowse) && (elem->GetNumChilds() > 0)) {
      fBrowsable.SetWorkingPath(path);
      return GetCurrentWorkingDirectory();
   }

   return ""s;
}

/////////////////////////////////////////////////////////////////////////////////
/// Show or update RBrowser in web window
/// If web window already started - just refresh it like "reload" button does
/// If no web window exists or \param always_start_new_browser configured, starts new window

void RBrowser::Show(const RWebDisplayArgs &args, bool always_start_new_browser)
{
   if (!fWebWindow->NumConnections() || always_start_new_browser) {
      fWebWindow->Show(args);
   } else {
      SendInitMsg(0);
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// Hide ROOT Browser

void RBrowser::Hide()
{
   if (!fWebWindow)
      return;

   fWebWindow->CloseConnections();
}


//////////////////////////////////////////////////////////////////////////////////////////////
/// Creates new widget

std::shared_ptr<RBrowserWidget> RBrowser::AddWidget(const std::string &kind)
{
   std::string name = kind + std::to_string(++fWidgetCnt);

   std::shared_ptr<RBrowserWidget> widget;

   if (kind == "editor")
      widget = std::make_shared<RBrowserEditorWidget>(name, true);
   else if (kind == "image")
      widget = std::make_shared<RBrowserEditorWidget>(name, false);
   else
      widget = RBrowserWidgetProvider::CreateWidget(kind, name);

   if (!widget) {
      R__LOG_ERROR(BrowserLog()) << "Fail to create widget of kind " << kind;
      return nullptr;
   }

   widget->Show("embed");
   fWidgets.emplace_back(widget);

   fActiveWidgetName = name;

   return widget;
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Returns active geometry viewer (if any)

std::shared_ptr<RBrowserWidget> RBrowser::FindWidget(const std::string &name) const
{
   auto iter = std::find_if(fWidgets.begin(), fWidgets.end(),
         [name](const std::shared_ptr<RBrowserWidget> &widget) { return name == widget->GetName(); });

   if (iter != fWidgets.end())
      return *iter;

   return nullptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Close and delete specified widget

void RBrowser::CloseTab(const std::string &name)
{
   auto iter = std::find_if(fWidgets.begin(), fWidgets.end(), [name](std::shared_ptr<RBrowserWidget> &widget) { return name == widget->GetName(); });
   if (iter != fWidgets.end())
      fWidgets.erase(iter);

   if (fActiveWidgetName == name)
      fActiveWidgetName.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Get content of history file

std::vector<std::string> RBrowser::GetRootHistory()
{
   std::vector<std::string> arr;

   std::string path = gSystem->UnixPathName(gSystem->HomeDirectory());
   path += "/.root_hist" ;
   std::ifstream infile(path);

   if (infile) {
      std::string line;
      while (std::getline(infile, line) && (arr.size() < 1000)) {
         if(!(std::find(arr.begin(), arr.end(), line) != arr.end())) {
            arr.emplace_back(line);
         }
      }
   }

   return arr;
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Get content of log file

std::vector<std::string> RBrowser::GetRootLogs()
{
   std::vector<std::string> arr;

   std::ostringstream pathtmp;
   pathtmp << gSystem->TempDirectory() << "/command." << gSystem->GetPid() << ".log";

   std::ifstream infile(pathtmp.str());
   if (infile) {
      std::string line;
      while (std::getline(infile, line) && (arr.size() < 10000)) {
         arr.emplace_back(line);
      }
   }

   return arr;
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Process client connect

void RBrowser::SendInitMsg(unsigned connid)
{
   std::vector<std::vector<std::string>> reply;

   reply.emplace_back(fBrowsable.GetWorkingPath()); // first element is current path

   for (auto &widget : fWidgets) {
      widget->ResetConn();
      reply.emplace_back(std::vector<std::string>({ widget->GetKind(), widget->GetUrl(), widget->GetName(), widget->GetTitle() }));
   }

   if (!fActiveWidgetName.empty())
      reply.emplace_back(std::vector<std::string>({ "active", fActiveWidgetName }));

   auto history = GetRootHistory();
   if (history.size() > 0) {
      history.insert(history.begin(), "history");
      reply.emplace_back(history);
   }

   auto logs = GetRootLogs();
   if (logs.size() > 0) {
      logs.insert(logs.begin(), "logs");
      reply.emplace_back(logs);
   }

   std::string msg = "INMSG:";
   msg.append(TBufferJSON::ToJSON(&reply, TBufferJSON::kNoSpaces).Data());

   fWebWindow->Send(connid, msg);
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Return the current directory of ROOT

std::string RBrowser::GetCurrentWorkingDirectory()
{
   return "WORKPATH:"s + TBufferJSON::ToJSON(&fBrowsable.GetWorkingPath()).Data();
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Create message which send to client to create new widget

std::string RBrowser::NewWidgetMsg(std::shared_ptr<RBrowserWidget> &widget)
{
   std::vector<std::string> arr = { widget->GetKind(), widget->GetUrl(), widget->GetName(), widget->GetTitle() };
   return "NEWWIDGET:"s + TBufferJSON::ToJSON(&arr, TBufferJSON::kNoSpaces).Data();
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Process received message from the client

void RBrowser::ProcessMsg(unsigned connid, const std::string &arg0)
{
   R__LOG_DEBUG(0, BrowserLog()) << "ProcessMsg  len " << arg0.length() << " substr(30) " << arg0.substr(0, 30);

   std::string kind, msg;
   auto pos = arg0.find(":");
   if (pos == std::string::npos) {
      kind = arg0;
   } else {
      kind = arg0.substr(0, pos);
      msg = arg0.substr(pos+1);
   }

   if (kind == "QUIT_ROOT") {

      fWebWindow->TerminateROOT();

   } else if (kind == "BRREQ") {
      // central place for processing browser requests
      auto json = ProcessBrowserRequest(msg);
      if (!json.empty()) fWebWindow->Send(connid, json);

   } else if (kind == "DBLCLK") {

      std::string reply;

      auto arr = TBufferJSON::FromJSON<std::vector<std::string>>(msg);
      if (arr && (arr->size() > 2))
         reply = ProcessDblClick(*arr);

      if (!reply.empty())
         fWebWindow->Send(connid, reply);

   } else if (kind == "WIDGET_SELECTED") {
      fActiveWidgetName = msg;
      auto widget = GetActiveWidget();
      if (widget) {
         auto reply = widget->SendWidgetContent();
         if (!reply.empty()) fWebWindow->Send(connid, reply);
      }
   } else if (kind == "CLOSE_TAB") {
      CloseTab(msg);
   } else if (kind == "GETWORKPATH") {
      fWebWindow->Send(connid, GetCurrentWorkingDirectory());
   } else if (kind == "CHPATH") {
      auto path = TBufferJSON::FromJSON<Browsable::RElementPath_t>(msg);
      if (path) fBrowsable.SetWorkingPath(*path);
      fWebWindow->Send(connid, GetCurrentWorkingDirectory());
   } else if (kind == "CMD") {
      std::string sPrompt = "root []";
      std::ostringstream pathtmp;
      pathtmp << gSystem->TempDirectory() << "/command." << gSystem->GetPid() << ".log";
      TApplication *app = gROOT->GetApplication();
      if (app->InheritsFrom("TRint")) {
         sPrompt = ((TRint*)gROOT->GetApplication())->GetPrompt();
         Gl_histadd((char *)msg.c_str());
      }

      std::ofstream ofs(pathtmp.str(), std::ofstream::out | std::ofstream::app);
      ofs << sPrompt << msg;
      ofs.close();

      gSystem->RedirectOutput(pathtmp.str().c_str(), "a");
      gROOT->ProcessLine(msg.c_str());
      gSystem->RedirectOutput(0);
   } else if (kind == "GETHISTORY") {

      auto history = GetRootHistory();

      fWebWindow->Send(connid, "HISTORY:"s + TBufferJSON::ToJSON(&history, TBufferJSON::kNoSpaces).Data());
   } else if (kind == "GETLOGS") {

      auto logs = GetRootLogs();
      fWebWindow->Send(connid, "LOGS:"s + TBufferJSON::ToJSON(&logs, TBufferJSON::kNoSpaces).Data());

   } else if (kind == "FILEDIALOG") {
      RFileDialog::Embedded(fWebWindow, arg0);
   } else if (kind == "SYNCEDITOR") {
      auto arr = TBufferJSON::FromJSON<std::vector<std::string>>(msg);
      if (arr && (arr->size() > 4)) {
         auto editor = std::dynamic_pointer_cast<RBrowserEditorWidget>(FindWidget(arr->at(0)));
         if (editor) {
            editor->fFirstSend = true;
            editor->fTitle = arr->at(1);
            editor->fFileName = arr->at(2);
            if (!arr->at(3).empty()) editor->fContent = arr->at(4);
            if ((arr->size() == 6) && (arr->at(5) == "SAVE"))
               ProcessSaveFile(editor->fFileName, editor->fContent);
            if ((arr->size() == 6) && (arr->at(5) == "RUN")) {
               ProcessSaveFile(editor->fFileName, editor->fContent);
               ProcessRunMacro(editor->fFileName);
            }

         }
      }
   } else if (kind == "NEWWIDGET") {
      auto widget = AddWidget(msg);
      if (widget)
         fWebWindow->Send(connid, NewWidgetMsg(widget));
   }
}
