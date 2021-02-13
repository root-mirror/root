/// \file gui_handler.cxx
// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2017-06-29
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#if !defined(_MSC_VER)
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wshadow"
#endif

#include "gui_handler.h"
#include "simple_app.h"

#include <sstream>
#include <string>

#include "include/base/cef_bind.h"
#include "include/cef_app.h"
#include "include/cef_version.h"
#include "include/views/cef_browser_view.h"
#include "include/views/cef_window.h"
#include "include/wrapper/cef_closure_task.h"
#include "include/wrapper/cef_helpers.h"
#include "include/cef_parser.h"
#include "include/wrapper/cef_stream_resource_handler.h"

#include "TEnv.h"
#include "TUrl.h"
#include "THttpServer.h"
#include "THttpCallArg.h"
#include "TSystem.h"
#include "TBase64.h"
#include <ROOT/RLogger.hxx>


ROOT::Experimental::RLogChannel &CefWebDisplayLog()
{
   static ROOT::Experimental::RLogChannel sChannel("ROOT.CefWebDisplay");
   return sChannel;
}


GuiHandler::GuiHandler(bool use_views) : fUseViews(use_views), is_closing_(false)
{
   fConsole = gEnv->GetValue("WebGui.Console", (int)0);

   // see https://bitbucket.org/chromiumembedded/cef-project/src/master/examples/resource_manager/?at=master for demo
   // one probably can avoid to use scheme handler and just redirect requests
   fResourceManager = new CefResourceManager();
}

void GuiHandler::OnTitleChange(CefRefPtr<CefBrowser> browser, const CefString &title)
{
   CEF_REQUIRE_UI_THREAD();

   if (fUseViews) {
      // Set the title of the window using the Views framework.
      CefRefPtr<CefBrowserView> browser_view = CefBrowserView::GetForBrowser(browser);
      if (browser_view) {
         CefRefPtr<CefWindow> window = browser_view->GetWindow();
         if (window) window->SetTitle(title);
      }
   } else {
      // Set the title of the window using platform APIs.
      PlatformTitleChange(browser, title);
   }
}

void GuiHandler::OnAfterCreated(CefRefPtr<CefBrowser> browser)
{
   CEF_REQUIRE_UI_THREAD();

   // Add to the list of existing browsers.
   fBrowserList.emplace_back(browser);
}

bool GuiHandler::DoClose(CefRefPtr<CefBrowser> browser)
{
   CEF_REQUIRE_UI_THREAD();

   // Closing the main window requires special handling. See the DoClose()
   // documentation in the CEF header for a detailed description of this
   // process.
   if (fBrowserList.size() == 1) {
      // Set a flag to indicate that the window close should be allowed.
      is_closing_ = true;
   }

   // Allow the close. For windowed browsers this will result in the OS close
   // event being sent.
   return false;
}

void GuiHandler::OnBeforeClose(CefRefPtr<CefBrowser> browser)
{
   CEF_REQUIRE_UI_THREAD();

   // Remove from the list of existing browsers.
   auto bit = fBrowserList.begin();
   for (; bit != fBrowserList.end(); ++bit) {
      if ((*bit)->IsSame(browser)) {
         fBrowserList.erase(bit);
         break;
      }
   }

   if (fBrowserList.empty()) {

      // All browser windows have closed. Quit the application message loop.

      CefQuitMessageLoop();
   }
}

// Returns a data: URI with the specified contents.
std::string GuiHandler::GetDataURI(const std::string& data, const std::string& mime_type)
{
    return "data:" + mime_type + ";base64," +
           CefURIEncode(CefBase64Encode(data.data(), data.size()), false)
            .ToString();
}


void GuiHandler::OnLoadError(CefRefPtr<CefBrowser> browser, CefRefPtr<CefFrame> frame, ErrorCode errorCode,
                              const CefString &errorText, const CefString &failedUrl)
{
   CEF_REQUIRE_UI_THREAD();

   // Don't display an error for downloaded files.
   if (errorCode == ERR_ABORTED)
      return;

   // Display a load error message.
   std::stringstream ss;
   ss << "<html><body bgcolor=\"white\">"
         "<h2>Failed to load URL "
      << failedUrl.ToString().substr(0,100) << " with error " << errorText.ToString() << " (" << errorCode
      << ").</h2></body></html>";
   // frame->LoadURL(GetDataURI(ss.str(), "text/html"));

   printf("Fail to load URL %s\n", failedUrl.ToString().substr(0,100).c_str());
}

void GuiHandler::CloseAllBrowsers(bool force_close)
{
   if (!CefCurrentlyOn(TID_UI)) {
      // Execute on the UI thread.
      CefPostTask(TID_UI, base::Bind(&GuiHandler::CloseAllBrowsers, this, force_close));
      return;
   }

   if (fBrowserList.empty())
      return;

   for (auto &br : fBrowserList)
      br->GetHost()->CloseBrowser(force_close);
}

bool GuiHandler::OnConsoleMessage(CefRefPtr<CefBrowser> browser,
                                  cef_log_severity_t level,
                                  const CefString &message, const CefString &source,
                                  int line)
{
   std::string src = source.ToString().substr(0,100);

   switch (level) {
   case LOGSEVERITY_WARNING:
      if (fConsole > -1)
         R__LOG_WARNING(CefWebDisplayLog()) << Form("CEF: %s:%d: %s", src.c_str(), line, message.ToString().c_str());
      break;
   case LOGSEVERITY_ERROR:
      if (fConsole > -2)
         R__LOG_ERROR(CefWebDisplayLog()) << Form("CEF: %s:%d: %s", src.c_str(), line, message.ToString().c_str());
      break;
   default:
      if (fConsole > 0)
         R__LOG_DEBUG(0, CefWebDisplayLog()) << Form("CEF: %s:%d: %s", src.c_str(), line, message.ToString().c_str());
      break;
   }

   return true;
}

cef_return_value_t GuiHandler::OnBeforeResourceLoad(
    CefRefPtr<CefBrowser> browser,
    CefRefPtr<CefFrame> frame,
    CefRefPtr<CefRequest> request,
    CefRefPtr<CefRequestCallback> callback) {
  CEF_REQUIRE_IO_THREAD();

  // std::string url = request->GetURL().ToString();
  // printf("OnBeforeResourceLoad url %s\n", url.c_str());

  return fResourceManager->OnBeforeResourceLoad(browser, frame, request,
                                                 callback);
}



class TCefHttpCallArg : public THttpCallArg {
protected:

   CefRefPtr<CefCallback> fCallBack{nullptr};

   void CheckWSPageContent(THttpWSHandler *) override
   {
      std::string search = "JSROOT.connectWebWindow({";
      std::string replace = search + "platform:\"cef3\",socket_kind:\"longpoll\",";

      ReplaceAllinContent(search, replace, true);
   }

public:
   explicit TCefHttpCallArg() = default;

   void AssignCallback(CefRefPtr<CefCallback> cb) { fCallBack = cb; }

   // this is callback from HTTP server
   void HttpReplied() override
   {
      if (IsFile()) {
         // send file
         std::string file_content = THttpServer::ReadFileContent((const char *)GetContent());
         SetContent(std::move(file_content));
      }

      fCallBack->Continue(); // we have result and can continue with processing
   }
};


class TGuiResourceHandler : public CefResourceHandler {
public:
   // QWebEngineUrlRequestJob *fRequest;

   THttpServer *fServer{nullptr};
   std::shared_ptr<TCefHttpCallArg> fArg;

   int fTransferOffset{0};

   explicit TGuiResourceHandler(THttpServer *serv, bool dummy = false)
   {
      fServer = serv;

      if (!dummy)
         fArg = std::make_shared<TCefHttpCallArg>();
   }

   virtual ~TGuiResourceHandler() {}

   void Cancel() OVERRIDE { CEF_REQUIRE_IO_THREAD(); }

   bool ProcessRequest(CefRefPtr<CefRequest> request, CefRefPtr<CefCallback> callback) OVERRIDE
   {
      CEF_REQUIRE_IO_THREAD();

      if (fArg) {
         fArg->AssignCallback(callback);
         fServer->SubmitHttp(fArg);
      } else {
         callback->Continue();
      }

      return true;
   }

   void GetResponseHeaders(CefRefPtr<CefResponse> response, int64 &response_length, CefString &redirectUrl) OVERRIDE
   {
      CEF_REQUIRE_IO_THREAD();

      if (!fArg || fArg->Is404()) {
         response->SetMimeType("text/html");
         response->SetStatus(404);
         response_length = 0;
      } else {
         response->SetMimeType(fArg->GetContentType());
         response->SetStatus(200);
         response_length = fArg->GetContentLength();

         if (fArg->NumHeader() > 0) {
            // printf("******* Response with extra headers\n");
            CefResponse::HeaderMap headers;
            for (Int_t n = 0; n < fArg->NumHeader(); ++n) {
               TString name = fArg->GetHeaderName(n);
               TString value = fArg->GetHeader(name.Data());
               headers.emplace(CefString(name.Data()), CefString(value.Data()));
               // printf("   header %s %s\n", name.Data(), value.Data());
            }
            response->SetHeaderMap(headers);
         }
//         if (strstr(fArg->GetQuery(),"connection="))
//            printf("Reply %s %s %s  len: %d %s\n", fArg->GetPathName(), fArg->GetFileName(), fArg->GetQuery(),
//                  fArg->GetContentLength(), (const char *) fArg->GetContent() );
      }
      // DCHECK(!fArg->Is404());
   }

   bool ReadResponse(void *data_out, int bytes_to_read, int &bytes_read, CefRefPtr<CefCallback> callback) OVERRIDE
   {
      CEF_REQUIRE_IO_THREAD();

      if (!fArg) return false;

      bytes_read = 0;

      if (fTransferOffset < fArg->GetContentLength()) {
         char *data_ = (char *)fArg->GetContent();
         // Copy the next block of data into the buffer.
         int transfer_size = fArg->GetContentLength() - fTransferOffset;
         if (transfer_size > bytes_to_read)
            transfer_size = bytes_to_read;
         memcpy(data_out, data_ + fTransferOffset, transfer_size);
         fTransferOffset += transfer_size;

         bytes_read = transfer_size;
      }

      // if content fully copied - can release reference, object will be cleaned up
      if (fTransferOffset >= fArg->GetContentLength())
         fArg.reset();

      return bytes_read > 0;
   }

   IMPLEMENT_REFCOUNTING(TGuiResourceHandler);
   DISALLOW_COPY_AND_ASSIGN(TGuiResourceHandler);
};


CefRefPtr<CefResourceHandler> GuiHandler::GetResourceHandler(
    CefRefPtr<CefBrowser> browser,
    CefRefPtr<CefFrame> frame,
    CefRefPtr<CefRequest> request) {
  CEF_REQUIRE_IO_THREAD();

  std::string addr = request->GetURL().ToString();
  std::string prefix = "http://rootserver.local";

  if (addr.compare(0, prefix.length(), prefix) != 0)
     return fResourceManager->GetResourceHandler(browser, frame, request);

  int indx = std::stoi(addr.substr(prefix.length(), addr.find("/", prefix.length()) - prefix.length()));

  if ((indx < 0) || (indx >= (int) fServers.size()) || !fServers[indx]) {
     R__LOG_ERROR(CefWebDisplayLog()) << "No THttpServer with index " << indx;
     return nullptr;
  }

  THttpServer *serv = fServers[indx];
  if (serv->IsZombie()) {
     fServers[indx] = nullptr;
     R__LOG_ERROR(CefWebDisplayLog()) << "THttpServer with index " << indx << " is zombie now";
     return nullptr;
  }

  TUrl url(addr.c_str());

  const char *inp_path = url.GetFile();

  TString inp_query = url.GetOptions();

  TString fname;

  if (serv->IsFileRequested(inp_path, fname)) {
     // process file - no need for special requests handling

     // when file not exists - return nullptr
     if (gSystem->AccessPathName(fname.Data()))
        return new TGuiResourceHandler(serv, true);

     const char *mime = THttpServer::GetMimeType(fname.Data());

     CefRefPtr<CefStreamReader> stream = CefStreamReader::CreateForFile(fname.Data());

     // Constructor for HTTP status code 200 and no custom response headers.
     // There’s also a version of the constructor for custom status code and response headers.
     return new CefStreamResourceHandler(mime, stream);
  }

  std::string inp_method = request->GetMethod().ToString();

  // printf("REQUEST METHOD %s\n", inp_method.c_str());

  TGuiResourceHandler *handler = new TGuiResourceHandler(serv);
  handler->fArg->SetMethod(inp_method.c_str());
  handler->fArg->SetPathAndFileName(inp_path);
  handler->fArg->SetTopName("webgui");

  if (inp_method == "POST") {

     CefRefPtr< CefPostData > post_data = request->GetPostData();

     if (!post_data) {
        R__LOG_ERROR(CefWebDisplayLog()) << "FATAL - NO POST DATA in CEF HANDLER!!!";
        exit(1);
     } else {
        CefPostData::ElementVector elements;
        post_data->GetElements(elements);
        size_t sz = 0, off = 0;
        for (unsigned n = 0; n < elements.size(); ++n)
           sz += elements[n]->GetBytesCount();
        std::string data;
        data.resize(sz);

        for (unsigned n = 0; n < elements.size(); ++n) {
           sz = elements[n]->GetBytes(elements[n]->GetBytesCount(), (char *)data.data() + off);
           off += sz;
        }
        handler->fArg->SetPostData(std::move(data));
     }
  } else if (inp_query.Index("&post=") != kNPOS) {
     Int_t pos = inp_query.Index("&post=");
     TString buf = TBase64::Decode(inp_query.Data() + pos + 6);
     handler->fArg->SetPostData(std::string(buf.Data()));
     inp_query.Resize(pos);
  }

  handler->fArg->SetQuery(inp_query.Data());

  // just return handler
  return handler;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Generate URL for batch page
/// Uses file:/// prefix to let access JSROOT scripts placed on file system
/// Register provider for that page in resource manager

std::string GuiHandler::AddBatchPage(const std::string &cont)
{
   std::string url = "file:///batch_page";
   url.append(std::to_string(fBatchPageCount++));
   url.append(".html");

   fResourceManager->AddContentProvider(url, cont, "text/html", 0, std::string());
   return url;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Generate URL for RWebWindow page
/// Register server instance and assign it with the index
/// Produced URL only works inside CEF and does not represent real HTTP address

std::string GuiHandler::MakePageUrl(THttpServer *serv, const std::string &addr)
{
   unsigned indx = 0;
   while ((indx < fServers.size()) && (fServers[indx] != serv)) indx++;
   if (indx >= fServers.size()) {
      fServers.emplace_back(serv);
      indx = fServers.size() - 1;
   }

   return std::string("http://rootserver.local") + std::to_string(indx) + addr;
}

