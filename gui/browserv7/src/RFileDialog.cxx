// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2019-10-31
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RFileDialog.hxx>

#include <ROOT/Browsable/RGroup.hxx>
#include <ROOT/Browsable/RSysFile.hxx>

#include <ROOT/RLogger.hxx>

#include "TBufferJSON.h"

#include <sstream>
#include <iostream>
#include <algorithm>
#include <memory>
#include <mutex>
#include <thread>
#include <fstream>

using namespace ROOT::Experimental;
using namespace std::string_literals;

/** \class ROOT::Experimental::RFileDialog
\ingroup rbrowser

web-based FileDialog.
*/

//////////////////////////////////////////////////////////////////////////////////////////////
/// constructor
/// When title not specified, default will be used

RFileDialog::RFileDialog(EDialogTypes kind, const std::string &title, const std::string &fname)
{
   fKind = kind;
   fTitle = title;

   if (fTitle.empty())
      switch (fKind) {
         case kOpenFile: fTitle = "Open file"; break;
         case kSaveAs: fTitle = "Save as file"; break;
         case kNewFile: fTitle = "New file"; break;
      }

   fSelect = fname;

   auto separ = fSelect.rfind("/");
   if (separ == std::string::npos)
      separ = fSelect.rfind("\\");

   std::string workdir;

   if (separ != std::string::npos) {
      workdir = fSelect.substr(0, separ);
      fSelect = fSelect.substr(separ+1);
   }

   auto comp = std::make_shared<Browsable::RGroup>("top", "Top file dialog element");
   auto workpath = Browsable::RSysFile::ProvideTopEntries(comp, workdir);
   fBrowsable.SetTopElement(comp);
   fBrowsable.SetWorkingPath(workpath);

   fWebWindow = RWebWindow::Create();

   // when dialog used in standalone mode, ui5 panel will be loaded
   fWebWindow->SetPanelName("rootui5.browser.view.FileDialog");

   // this is call-back, invoked when message received via websocket
   fWebWindow->SetCallBacks([this](unsigned connid) { SendInitMsg(connid); },
                            [this](unsigned connid, const std::string &arg) { ProcessMsg(connid, arg); },
                            [this](unsigned) { InvokeCallBack(); });
   fWebWindow->SetGeometry(800, 600); // configure predefined window geometry
   fWebWindow->SetConnLimit(1); // the only connection is allowed
   fWebWindow->SetMaxQueueLength(30); // number of allowed entries in the window queue
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// destructor

RFileDialog::~RFileDialog()
{
   InvokeCallBack(); // invoke callback if not yet performed
   R__LOG_DEBUG(0, BrowserLog()) << "RFileDialog destructor";
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Assign callback. If file was already selected, immediately call it

void RFileDialog::SetCallback(RFileDialogCallback_t callback)
{
   fCallback = callback;
   if (fDidSelect)
      InvokeCallBack();
}

/////////////////////////////////////////////////////////////////////////////////
/// Show or update RFileDialog in web window
/// If web window already started - just refresh it like "reload" button does
/// Reset result of file selection (if any)

void RFileDialog::Show(const RWebDisplayArgs &args)
{
   fDidSelect = false;

   if (fWebWindow->NumConnections() == 0) {
      RWebWindow::ShowWindow(fWebWindow, args);
   } else {
      SendInitMsg(0);
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// Hide ROOT Browser

void RFileDialog::Hide()
{
   fWebWindow->CloseConnections();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns dialog type as string
/// String value used for configuring JS-side

std::string RFileDialog::TypeAsString(EDialogTypes kind)
{
   switch(kind) {
      case kOpenFile: return "OpenFile"s;
      case kSaveAs: return "SaveAs"s;
      case kNewFile: return "NewFile"s;
   }

   return ""s;
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Configure selected filter
/// Has to be one of the string from NameFilters entry

void RFileDialog::SetSelectedFilter(const std::string &name)
{
   fSelectedFilter = name;
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Returns selected filter
/// Can differ from specified value - if it does not match to existing entry in NameFilters

std::string RFileDialog::GetSelectedFilter() const
{
   if (fNameFilters.size() == 0)
      return fSelectedFilter;

   std::string lastname, allname;

   for (auto &entry : fNameFilters) {
      auto pp = entry.find(" (");
      if (pp == std::string::npos) continue;
      auto name = entry.substr(0, pp);

      if (name == fSelectedFilter)
         return name;

      if (allname.empty() && GetRegexp(name).empty())
         allname = name;

      lastname = name;
   }

   if (!allname.empty()) return allname;
   if (!lastname.empty()) return lastname;

   return ""s;
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Returns regexp for selected filter
/// String should have form "Filter name (*.ext1 *.ext2 ...)

std::string RFileDialog::GetRegexp(const std::string &fname) const
{
   if (!fname.empty())
      for (auto &entry : fNameFilters) {
         if (entry.compare(0, fname.length(), fname) == 0) {
            auto pp = entry.find("(", fname.length());

            std::string res;

            while (pp != std::string::npos) {
               pp = entry.find("*.", pp);
               if (pp == std::string::npos) break;

               auto pp2 = entry.find_first_of(" )", pp+2);
               if (pp2 == std::string::npos) break;

               if (res.empty()) res = "^(.*\\.(";
                           else res.append("|");

               res.append(entry.substr(pp+2, pp2 - pp - 2));

               pp = pp2;
            }

            if (!res.empty()) res.append(")$)");

            return res;

         }
      }

   return ""s;
}


//////////////////////////////////////////////////////////////////////////////////////////////
/// Sends initial message to the client

void RFileDialog::SendInitMsg(unsigned connid)
{
   auto filter = GetSelectedFilter();
   RBrowserRequest req;
   req.sort = "alphabetical";
   req.regex = GetRegexp(filter);

   auto jtitle = TBufferJSON::ToJSON(&fTitle);
   auto jpath = TBufferJSON::ToJSON(&fBrowsable.GetWorkingPath());
   auto jfname = TBufferJSON::ToJSON(&fSelect);
   auto jfilters = TBufferJSON::ToJSON(&fNameFilters);
   auto jfilter = TBufferJSON::ToJSON(&filter);

   fWebWindow->Send(connid, "INMSG:{\"kind\" : \""s + TypeAsString(fKind) + "\", "s +
                                   "\"title\" : "s + jtitle.Data() + ","s +
                                   "\"path\" : "s + jpath.Data() + ","s +
                                   "\"filter\" : "s + jfilter.Data() + ","s +
                                   "\"filters\" : "s + jfilters.Data() + ","s +
                                   "\"fname\" : "s + jfname.Data() + ","s +
                                   "\"brepl\" : "s + fBrowsable.ProcessRequest(req) + " }"s);
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Sends new data after change current directory

void RFileDialog::SendChPathMsg(unsigned connid)
{
   RBrowserRequest req;
   req.sort = "alphabetical";
   req.regex = GetRegexp(GetSelectedFilter());

   auto jpath = TBufferJSON::ToJSON(&fBrowsable.GetWorkingPath());

   fWebWindow->Send(connid, "CHMSG:{\"path\" : "s + jpath.Data() +
                                 ", \"brepl\" : "s + fBrowsable.ProcessRequest(req) + " }"s);
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Process received data from client

void RFileDialog::ProcessMsg(unsigned connid, const std::string &arg)
{
   if (arg.compare(0, 7, "CHPATH:") == 0) {
      auto path = TBufferJSON::FromJSON<Browsable::RElementPath_t>(arg.substr(7));
      if (path) fBrowsable.SetWorkingPath(*path);

      SendChPathMsg(connid);

   } else if (arg.compare(0, 6, "CHEXT:") == 0) {

      SetSelectedFilter(arg.substr(6));

      SendChPathMsg(connid);

   } else if (arg.compare(0, 10, "DLGSELECT:") == 0) {
      // selected file name, if file exists - send request for confirmation

      auto path = TBufferJSON::FromJSON<Browsable::RElementPath_t>(arg.substr(10));

      if (!path) {
         R__LOG_ERROR(BrowserLog()) << "Fail to decode JSON " << arg.substr(10);
         return;
      }

      // check if element exists
      auto elem = fBrowsable.GetElementFromTop(*path);
      if (elem)
         fSelect = elem->GetContent("filename");
      else
         fSelect.clear();

      bool need_confirm = false;

      if ((GetType() == kSaveAs) || (GetType() == kNewFile)) {
         if (elem) {
            need_confirm = true;
         } else {
            std::string fname = path->back();
            path->pop_back();
            auto direlem = fBrowsable.GetElementFromTop(*path);
            if (direlem)
               fSelect = direlem->GetContent("filename") + "/"s + fname;
         }
      }

      if (need_confirm) {
         fWebWindow->Send(connid, "NEED_CONFIRM"s); // sending request for confirmation
      } else {
         fWebWindow->Send(connid, "SELECT_CONFIRMED:"s + fSelect); // sending select confirmation with fully qualified file name
         fDidSelect = true;
      }
   } else if (arg == "DLGNOSELECT") {
      fSelect.clear();
      fDidSelect = true;
      fWebWindow->Send(connid, "NOSELECT_CONFIRMED"s); // sending confirmation of NOSELECT
   } else if (arg == "DLG_CONFIRM_SELECT") {
      fDidSelect = true;
      fWebWindow->Send(connid, "SELECT_CONFIRMED:"s + fSelect);
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Invoke specified callback

void RFileDialog::InvokeCallBack()
{
   if (fCallback) {
      auto func = fCallback;
      // reset callback to release associated with lambda resources
      // reset before invoking callback to avoid multiple calls
      fCallback = nullptr;
      func(fSelect);
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Start specified dialog type

std::string RFileDialog::Dialog(EDialogTypes kind, const std::string &title, const std::string &fname)
{
   RFileDialog dlg(kind, title, fname);

   dlg.Show();

   dlg.fWebWindow->WaitForTimed([&](double) {
      if (dlg.fDidSelect) return 1;

      return 0; // continue waiting
   });

   return dlg.fSelect;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Start OpenFile dialog.
/// Blocks until file name is selected or Cancel button is pressed
/// Returns selected file name (or empty string)

std::string RFileDialog::OpenFile(const std::string &title, const std::string &fname)
{
   return Dialog(kOpenFile, title, fname);
}

/////////////////////////////////////////////////////////////////////////////////////
/// Start SaveAs dialog.
/// Blocks until file name is selected or Cancel button is pressed
/// Returns selected file name (or empty string)

std::string RFileDialog::SaveAs(const std::string &title, const std::string &fname)
{
   return Dialog(kSaveAs, title, fname);
}

/////////////////////////////////////////////////////////////////////////////////////
/// Start NewFile dialog.
/// Blocks until file name is selected or Cancel button is pressed
/// Returns selected file name (or empty string)

std::string RFileDialog::NewFile(const std::string &title, const std::string &fname)
{
   return Dialog(kNewFile, title, fname);
}

/////////////////////////////////////////////////////////////////////////////////////
/// Create dialog instance to use as embedded dialog inside other widget
/// Embedded dialog started on the client side where FileDialogController.SaveAs() method called
/// Such method immediately send message with "FILEDIALOG:" prefix
/// On the server side widget should detect such message and call RFileDialog::Embedded()
/// providing received string as second argument.
/// Returned instance of shared_ptr<RFileDialog> may be used to assign callback when file is selected

std::shared_ptr<RFileDialog> RFileDialog::Embedded(const std::shared_ptr<RWebWindow> &window, const std::string &args)
{
   if (args.compare(0, 11, "FILEDIALOG:") != 0)
      return nullptr;

   auto arr = TBufferJSON::FromJSON<std::vector<std::string>>(args.substr(11));

   if (!arr || (arr->size() < 3)) {
      R__LOG_ERROR(BrowserLog()) << "Embedded FileDialog failure - argument should have at least three strings" << args.substr(11);
      return nullptr;
   }

   auto kind = kSaveAs;

   if (TypeAsString(kOpenFile) == arr->at(0))
      kind = kOpenFile;
   else if (TypeAsString(kNewFile) == arr->at(0))
      kind = kNewFile;

   auto dialog = std::make_shared<RFileDialog>(kind, "", arr->at(1));

   auto chid = std::stoi(arr->at(2));

   if (arr->size() > 4) {
      dialog->SetSelectedFilter(arr->at(3));
      arr->erase(arr->begin(), arr->begin() + 4); // erase 4 elements, keep only list of filters
      dialog->SetNameFilters(*arr);
   }

   dialog->Show({window, chid});

   // use callback to release pointer, actually not needed but just to avoid compiler warning
   dialog->SetCallback([dialog](const std::string &) mutable { dialog.reset(); });

   return dialog;
}
