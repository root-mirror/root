// Author: Sergey Linev <s.linev@gsi.de>
// Date: 2018-10-17
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RWebDisplayHandle.hxx>

#include <ROOT/RMakeUnique.hxx>
#include <ROOT/RLogger.hxx>

#include "RConfigure.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TString.h"
#include "TObjArray.h"
#include "THttpServer.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TBase64.h"

#include <regex>
#include <fstream>

#ifdef _MSC_VER
#include <process.h>
#else
#include <unistd.h>
#include <stdlib.h>
#include <signal.h>
#include <spawn.h>
#endif

using namespace ROOT::Experimental;
using namespace std::string_literals;

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Static holder of registered creators of web displays

std::map<std::string, std::unique_ptr<RWebDisplayHandle::Creator>> &RWebDisplayHandle::GetMap()
{
   static std::map<std::string, std::unique_ptr<RWebDisplayHandle::Creator>> sMap;
   return sMap;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Search for specific browser creator
/// If not found, try to add one
/// \param name - creator name like ChromeCreator
/// \param libname - shared library name where creator could be provided

std::unique_ptr<RWebDisplayHandle::Creator> &RWebDisplayHandle::FindCreator(const std::string &name, const std::string &libname)
{
   auto &m = GetMap();
   auto search = m.find(name);
   if (search == m.end()) {

      if (libname == "ChromeCreator") {
         m.emplace(name, std::make_unique<ChromeCreator>());
      } else if (libname == "FirefoxCreator") {
         m.emplace(name, std::make_unique<FirefoxCreator>());
      } else if (libname == "BrowserCreator") {
         m.emplace(name, std::make_unique<BrowserCreator>(false));
      } else if (!libname.empty()) {
         gSystem->Load(libname.c_str());
      }

      search = m.find(name); // try again
   }

   if (search != m.end())
      return search->second;

   static std::unique_ptr<RWebDisplayHandle::Creator> dummy;
   return dummy;
}

namespace ROOT {
namespace Experimental {

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Specialized handle to hold information about running browser process
/// Used to correctly cleanup all processes and temporary directories

class RWebBrowserHandle : public RWebDisplayHandle {

#ifdef _MSC_VER
   typedef int browser_process_id;
#else
   typedef pid_t browser_process_id;
#endif
   std::string fTmpDir;         ///< temporary directory to delete at the end
   bool fHasPid{false};
   browser_process_id fPid;

public:
   RWebBrowserHandle(const std::string &url, const std::string &tmpdir, const std::string &dump) : RWebDisplayHandle(url), fTmpDir(tmpdir)
   {
      SetContent(dump);
   }

   RWebBrowserHandle(const std::string &url, const std::string &tmpdir, browser_process_id pid)
      : RWebDisplayHandle(url), fTmpDir(tmpdir), fHasPid(true), fPid(pid)
   {
   }

   virtual ~RWebBrowserHandle()
   {
#ifdef _MSC_VER
      if (fHasPid)
         gSystem->Exec(("taskkill /F /PID "s + std::to_string(fPid) + " >NUL 2>NUL").c_str());
      std::string rmdir = "rmdir /S /Q ";
#else
      if (fHasPid)
         kill(fPid, SIGKILL);
      std::string rmdir = "rm -rf ";
#endif
      if (!fTmpDir.empty())
         gSystem->Exec((rmdir + fTmpDir).c_str());
   }

};

} // namespace Experimental
} // namespace ROOT

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Class to handle starting of web-browsers like Chrome or Firefox

RWebDisplayHandle::BrowserCreator::BrowserCreator(bool custom, const std::string &exec)
{
   if (custom) return;

   if (!exec.empty()) {
      if (exec.find("$url") == std::string::npos) {
         fProg = exec;
#ifdef _MSC_VER
         fExec = exec + " $url";
#else
         fExec = exec + " $url &";
#endif
      } else {
         fExec = exec;
         auto pos = exec.find(" ");
         if (pos != std::string::npos)
            fProg = exec.substr(0, pos);
      }
   } else if (gSystem->InheritsFrom("TMacOSXSystem")) {
      fExec = "open \'$url\'";
   } else if (gSystem->InheritsFrom("TWinNTSystem")) {
      fExec = "start $url";
   } else {
      fExec = "xdg-open \'$url\' &";
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Check if browser executable exists and can be used

void RWebDisplayHandle::BrowserCreator::TestProg(const std::string &nexttry, bool check_std_paths)
{
   if (nexttry.empty() || !fProg.empty())
      return;

   if (!gSystem->AccessPathName(nexttry.c_str(), kExecutePermission)) {
#ifdef R__MACOSX
      fProg = std::regex_replace(nexttry, std::regex("%20"), " ");
#else
      fProg = nexttry;
#endif
      return;
   }

   if (!check_std_paths)
      return;

#ifdef _MSC_VER
   std::string ProgramFiles = gSystem->Getenv("ProgramFiles");
   auto pos = ProgramFiles.find(" (x86)");
   if (pos != std::string::npos)
      ProgramFiles.erase(pos, 6);
   std::string ProgramFilesx86 = gSystem->Getenv("ProgramFiles(x86)");

   if (!ProgramFiles.empty())
      TestProg(ProgramFiles + nexttry, false);
   if (!ProgramFilesx86.empty())
      TestProg(ProgramFilesx86 + nexttry, false);
#endif
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Display given URL in web browser

std::unique_ptr<RWebDisplayHandle>
RWebDisplayHandle::BrowserCreator::Display(const RWebDisplayArgs &args)
{
   std::string url = args.GetFullUrl();
   if (url.empty())
      return nullptr;

   std::string exec;
   if (args.IsHeadless())
      exec = fBatchExec;
   else if (args.IsStandalone())
      exec = fExec;
   else
      exec = "$prog $url &";

   if (exec.empty())
      return nullptr;

   std::string swidth = std::to_string(args.GetWidth() > 0 ? args.GetWidth() : 800),
               sheight = std::to_string(args.GetHeight() > 0 ? args.GetHeight() : 600),
               sposx = std::to_string(args.GetX() >= 0 ? args.GetX() : 0),
               sposy = std::to_string(args.GetY() >= 0 ? args.GetY() : 0);

   ProcessGeometry(exec, args);

   std::string rmdir = MakeProfile(exec, args.IsHeadless());

   exec = std::regex_replace(exec, std::regex("\\$url"), url);
   exec = std::regex_replace(exec, std::regex("\\$width"), swidth);
   exec = std::regex_replace(exec, std::regex("\\$height"), sheight);
   exec = std::regex_replace(exec, std::regex("\\$posx"), sposx);
   exec = std::regex_replace(exec, std::regex("\\$posy"), sposy);

   if (exec.compare(0,5,"fork:") == 0) {
      if (fProg.empty()) {
         R__LOG_ERROR(WebGUILog()) << "Fork instruction without executable";
         return nullptr;
      }

      exec.erase(0, 5);

#ifndef _MSC_VER

      std::unique_ptr<TObjArray> fargs(TString(exec.c_str()).Tokenize(" "));
      if (!fargs || (fargs->GetLast()<=0)) {
         R__LOG_ERROR(WebGUILog()) << "Fork instruction is empty";
         return nullptr;
      }

      std::vector<char *> argv;
      argv.push_back((char *) fProg.c_str());
      for (Int_t n = 0; n <= fargs->GetLast(); ++n)
         argv.push_back((char *)fargs->At(n)->GetName());
      argv.push_back(nullptr);

      R__LOG_DEBUG(0, WebGUILog()) << "Show web window in browser with posix_spawn:\n" << fProg << " " << exec;

      pid_t pid;
      int status = posix_spawn(&pid, argv[0], nullptr, nullptr, argv.data(), nullptr);
      if (status != 0) {
         R__LOG_ERROR(WebGUILog()) << "Fail to launch " << argv[0];
         return nullptr;
      }

      // add processid and rm dir

      return std::make_unique<RWebBrowserHandle>(url, rmdir, pid);

#else

      if (fProg.empty()) {
         R__LOG_ERROR(WebGUILog()) << "No Web browser found";
         return nullptr;
      }

      // use UnixPathName to simplify handling of backslashes
      exec = "wmic process call create '"s + gSystem->UnixPathName(fProg.c_str()) + exec + "' | find \"ProcessId\" "s;
      std::string process_id = gSystem->GetFromPipe(exec.c_str());
      std::stringstream ss(process_id);
      std::string tmp;
      char c;
      int pid = 0;
      ss >> tmp >> c >> pid;

      if (pid <= 0) {
         R__LOG_ERROR(WebGUILog()) << "Fail to launch " << fProg;
         return nullptr;
      }

      // add processid and rm dir
      return std::make_unique<RWebBrowserHandle>(url, rmdir, pid);
#endif
   }

#ifdef _MSC_VER

   if (exec.rfind("&") == exec.length() - 1) {

      // if last symbol is &, use _spawn to detach execution
      exec.resize(exec.length() - 1);

      std::vector<char *> argv;
      std::string firstarg = fProg;
      auto slashpos = firstarg.find_last_of("/\\");
      if (slashpos != std::string::npos)
         firstarg.erase(0, slashpos + 1);
      argv.push_back((char *)firstarg.c_str());

      std::unique_ptr<TObjArray> fargs(TString(exec.c_str()).Tokenize(" "));
      for (Int_t n = 1; n <= fargs->GetLast(); ++n)
         argv.push_back((char *)fargs->At(n)->GetName());
      argv.push_back(nullptr);

      R__LOG_DEBUG(0, WebGUILog()) << "Showing web window in " << fProg << " with:\n" << exec;

      _spawnv(_P_NOWAIT, fProg.c_str(), argv.data());

      return std::make_unique<RWebBrowserHandle>(url, rmdir, ""s);
   }

   std::string prog = "\""s + gSystem->UnixPathName(fProg.c_str()) + "\""s;

#else

#ifdef R__MACOSX
   std::string prog = std::regex_replace(fProg, std::regex(" "), "\\ ");
#else
   std::string prog = fProg;
#endif

#endif

   exec = std::regex_replace(exec, std::regex("\\$prog"), prog);

   std::string redirect = args.GetRedirectOutput(), dump_content;

   if (!redirect.empty()) {
      auto p = exec.length();
      if (exec.rfind("&") == p-1) --p;
      exec.insert(p, " >"s + redirect + " "s);
   }

   R__LOG_DEBUG(0, WebGUILog()) << "Showing web window in browser with:\n" << exec;

   gSystem->Exec(exec.c_str());

   // read content of redirected output
   if (!redirect.empty()) {
      dump_content = THttpServer::ReadFileContent(redirect.c_str());

      gSystem->Unlink(redirect.c_str());
   }

   // add rmdir if required
   return std::make_unique<RWebBrowserHandle>(url, rmdir, dump_content);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor

RWebDisplayHandle::ChromeCreator::ChromeCreator() : BrowserCreator(true)
{
   TestProg(gEnv->GetValue("WebGui.Chrome", ""));

#ifdef _MSC_VER
   TestProg("\\Google\\Chrome\\Application\\chrome.exe", true);
#endif
#ifdef R__MACOSX
   TestProg("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome");
#endif
#ifdef R__LINUX
   TestProg("/usr/bin/chromium");
   TestProg("/usr/bin/chromium-browser");
   TestProg("/usr/bin/chrome-browser");
#endif

#ifdef _MSC_VER
   // fBatchExec = gEnv->GetValue("WebGui.ChromeBatch", "fork: --headless --disable-gpu $geometry $url");

   fBatchExec = gEnv->GetValue("WebGui.ChromeBatch", "$prog --headless $geometry $url");
   fExec = gEnv->GetValue("WebGui.ChromeInteractive", "$prog $geometry --no-first-run --app=$url &"); // & in windows mean usage of spawn
#else
   fBatchExec = gEnv->GetValue("WebGui.ChromeBatch", "$prog --headless --incognito $geometry $url");
   fExec = gEnv->GetValue("WebGui.ChromeInteractive", "$prog $geometry --no-first-run --incognito --app=\'$url\' &");
#endif
}


//////////////////////////////////////////////////////////////////////////////////////////////////
/// Replace $geometry placeholder with geometry settings
/// Also RWebDisplayArgs::GetExtraArgs() are appended

void RWebDisplayHandle::ChromeCreator::ProcessGeometry(std::string &exec, const RWebDisplayArgs &args)
{
   std::string geometry;
   if ((args.GetWidth() > 0) && (args.GetHeight() > 0))
      geometry = "--window-size="s + std::to_string(args.GetWidth())
                                   + (args.IsHeadless() ? "x"s : ","s)
                                   + std::to_string(args.GetHeight());

   if (((args.GetX() >= 0) || (args.GetY() >= 0)) && !args.IsHeadless()) {
      if (!geometry.empty()) geometry.append(" ");
      geometry.append("--window-position="s + std::to_string(args.GetX() >= 0 ? args.GetX() : 0) + ","s +
                                           std::to_string(args.GetY() >= 0 ? args.GetY() : 0));
   }

   if (!args.GetExtraArgs().empty()) {
      if (!geometry.empty()) geometry.append(" ");
      geometry.append(args.GetExtraArgs());
   }

   exec = std::regex_replace(exec, std::regex("\\$geometry"), geometry);
}


//////////////////////////////////////////////////////////////////////////////////////////////////
/// Handle profile argument

std::string RWebDisplayHandle::ChromeCreator::MakeProfile(std::string &exec, bool)
{
   std::string rmdir, profile_arg;

   if (exec.find("$profile") == std::string::npos)
      return rmdir;

   const char *chrome_profile = gEnv->GetValue("WebGui.ChromeProfile", "");
   if (chrome_profile && *chrome_profile) {
      profile_arg = chrome_profile;
   } else {
      gRandom->SetSeed(0);
      rmdir = profile_arg = std::string(gSystem->TempDirectory()) + "/root_chrome_profile_"s + std::to_string(gRandom->Integer(0x100000));
   }

   exec = std::regex_replace(exec, std::regex("\\$profile"), profile_arg);

   return rmdir;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor

RWebDisplayHandle::FirefoxCreator::FirefoxCreator() : BrowserCreator(true)
{
   TestProg(gEnv->GetValue("WebGui.Firefox", ""));

#ifdef _MSC_VER
   TestProg("\\Mozilla Firefox\\firefox.exe", true);
#endif
#ifdef R__MACOSX
   TestProg("/Applications/Firefox.app/Contents/MacOS/firefox");
#endif
#ifdef R__LINUX
   TestProg("/usr/bin/firefox");
#endif

#ifdef _MSC_VER
   // there is a problem when specifying the window size with wmic on windows:
   // It gives: Invalid format. Hint: <paramlist> = <param> [, <paramlist>].
   // fBatchExec = gEnv->GetValue("WebGui.FirefoxBatch", "fork: -headless -no-remote $profile $url");
   fExec = gEnv->GetValue("WebGui.FirefoxInteractive", "$prog -no-remote $profile $url &");
#else
   // fBatchExec = gEnv->GetValue("WebGui.FirefoxBatch", "fork:--headless --private-window --no-remote $profile $url");
   fExec = gEnv->GetValue("WebGui.FirefoxInteractive", "$prog --private-window \'$url\' &");
#endif
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Create Firefox profile to run independent browser window

std::string RWebDisplayHandle::FirefoxCreator::MakeProfile(std::string &exec, bool batch_mode)
{
   std::string rmdir, profile_arg;

   if (exec.find("$profile") == std::string::npos)
      return rmdir;

   const char *ff_profile = gEnv->GetValue("WebGui.FirefoxProfile", "");
   const char *ff_profilepath = gEnv->GetValue("WebGui.FirefoxProfilePath", "");
   Int_t ff_randomprofile = gEnv->GetValue("WebGui.FirefoxRandomProfile", (Int_t) 0);
   if (ff_profile && *ff_profile) {
      profile_arg = "-P "s + ff_profile;
   } else if (ff_profilepath && *ff_profilepath) {
      profile_arg = "-profile "s + ff_profilepath;
   } else if ((ff_randomprofile > 0) || (batch_mode && (ff_randomprofile >= 0))) {

      gRandom->SetSeed(0);
      std::string rnd_profile = "root_ff_profile_"s + std::to_string(gRandom->Integer(0x100000));
      std::string profile_dir = std::string(gSystem->TempDirectory()) + "/"s + rnd_profile;

      profile_arg = "-profile "s + profile_dir;

      if (gSystem->mkdir(profile_dir.c_str()) == 0) {
         rmdir = profile_dir;
      } else {
         R__LOG_ERROR(WebGUILog()) << "Cannot create Firefox profile directory " << profile_dir;
      }
   }

   exec = std::regex_replace(exec, std::regex("\\$profile"), profile_arg);

   return rmdir;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Create web display
/// \param args - defines where and how to display web window
/// Returns RWebDisplayHandle, which holds information of running browser application
/// Can be used fully independent from RWebWindow classes just to show any web page

std::unique_ptr<RWebDisplayHandle> RWebDisplayHandle::Display(const RWebDisplayArgs &args)
{
   std::unique_ptr<RWebDisplayHandle> handle;

   auto try_creator = [&](std::unique_ptr<Creator> &creator) {
      if (!creator || !creator->IsActive())
         return false;
      handle = creator->Display(args);
      return handle ? true : false;
   };

   if ((args.GetBrowserKind() == RWebDisplayArgs::kLocal) || (args.GetBrowserKind() == RWebDisplayArgs::kCEF)) {
      if (try_creator(FindCreator("cef", "libROOTCefDisplay")))
         return handle;
   }

   if ((args.GetBrowserKind() == RWebDisplayArgs::kLocal) || (args.GetBrowserKind() == RWebDisplayArgs::kQt5)) {
      if (try_creator(FindCreator("qt5", "libROOTQt5WebDisplay")))
         return handle;
   }

   if (args.IsLocalDisplay()) {
      R__LOG_ERROR(WebGUILog()) << "Neither Qt5 nor CEF libraries were found to provide local display";
      return handle;
   }

   if ((args.GetBrowserKind() == RWebDisplayArgs::kNative) || (args.GetBrowserKind() == RWebDisplayArgs::kChrome)) {
      if (try_creator(FindCreator("chrome", "ChromeCreator")))
         return handle;
   }

   if ((args.GetBrowserKind() == RWebDisplayArgs::kNative) || (args.GetBrowserKind() == RWebDisplayArgs::kFirefox)) {
      if (try_creator(FindCreator("firefox", "FirefoxCreator")))
         return handle;
   }

   if ((args.GetBrowserKind() == RWebDisplayArgs::kChrome) || (args.GetBrowserKind() == RWebDisplayArgs::kFirefox)) {
      // R__LOG_ERROR(WebGUILog()) << "Neither Chrome nor Firefox browser cannot be started to provide display";
      return handle;
   }

   if ((args.GetBrowserKind() == RWebDisplayArgs::kCustom)) {
      std::unique_ptr<Creator> creator = std::make_unique<BrowserCreator>(false, args.GetCustomExec());
      try_creator(creator);
   } else {
      try_creator(FindCreator("browser", "BrowserCreator"));
   }

   return handle;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Display provided url in configured web browser
/// \param url - specified URL address like https://root.cern
/// Browser can specified when starting `root --web=firefox`
/// Returns true when browser started
/// It is convenience method, equivalent to:
///  ~~~
///     RWebDisplayArgs args;
///     args.SetUrl(url);
///     args.SetStandalone(false);
///     auto handle = RWebDisplayHandle::Display(args);
/// ~~~

bool RWebDisplayHandle::DisplayUrl(const std::string &url)
{
   RWebDisplayArgs args;
   args.SetUrl(url);
   args.SetStandalone(false);

   auto handle = Display(args);

   return !!handle;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
/// Produce image file using JSON data as source
/// Invokes JSROOT drawing functionality in headless browser - Google Chrome

bool RWebDisplayHandle::ProduceImage(const std::string &fname, const std::string &json, int width, int height)
{
   if (json.empty())
      return false;

   std::string _fname = fname;
   std::transform(_fname.begin(), _fname.end(), _fname.begin(), ::tolower);

   auto EndsWith = [_fname](const std::string &suffix) {
      return (_fname.length() > suffix.length()) ? (0 == _fname.compare (_fname.length() - suffix.length(), suffix.length(), suffix)) : false;
   };

   if (EndsWith(".json")) {
      std::ofstream ofs(fname);
      ofs << json;
      return true;
   }

   const char *jsrootsys = gSystem->Getenv("JSROOTSYS");
   TString jsrootsysdflt;
   if (!jsrootsys) {
      jsrootsysdflt = TROOT::GetDataDir() + "/js";
      if (gSystem->ExpandPathName(jsrootsysdflt)) {
         R__LOG_ERROR(WebGUILog()) << "Fail to locate JSROOT " << jsrootsysdflt;
         return false;
      }
      jsrootsys = jsrootsysdflt.Data();
   }

   RWebDisplayArgs args; // set default browser kind, only Chrome or CEF or Qt5 can be used here
   if ((args.GetBrowserKind() != RWebDisplayArgs::kCEF) && (args.GetBrowserKind() != RWebDisplayArgs::kQt5))
      args.SetBrowserKind(RWebDisplayArgs::kChrome);

   std::string draw_kind;

   if (EndsWith(".pdf"))
      draw_kind = "draw"; // not a JSROOT drawing but Chrome capability to create PDF out of HTML page is used
   else if (EndsWith("shot.png"))
      draw_kind = (args.GetBrowserKind() == RWebDisplayArgs::kChrome) ? "draw" : "png";
   else if (EndsWith(".svg"))
      draw_kind = "svg";
   else if (EndsWith(".png"))
      draw_kind = "png";
   else if (EndsWith(".jpg") || EndsWith(".jpeg"))
      draw_kind = "jpeg";
   else if (EndsWith(".webp"))
      draw_kind = "webp";
   else
      return false;

   TString origin = TROOT::GetDataDir() + "/js/files/canv_batch.htm";
   if (gSystem->ExpandPathName(origin)) {
      R__LOG_ERROR(WebGUILog()) << "Fail to find " << origin;
      return false;
   }

   auto filecont = THttpServer::ReadFileContent(origin.Data());
   if (filecont.empty()) {
      R__LOG_ERROR(WebGUILog()) << "Fail to read content of " << origin;
      return false;
   }

   filecont = std::regex_replace(filecont, std::regex("\\$draw_width"), std::to_string(width));
   filecont = std::regex_replace(filecont, std::regex("\\$draw_height"), std::to_string(height));

   if (strstr(jsrootsys,"http://") || strstr(jsrootsys,"https://") || strstr(jsrootsys,"file://"))
      filecont = std::regex_replace(filecont, std::regex("\\$jsrootsys"), jsrootsys);
   else
      filecont = std::regex_replace(filecont, std::regex("\\$jsrootsys"), "file://"s + jsrootsys);

   filecont = std::regex_replace(filecont, std::regex("\\$draw_kind"), draw_kind);

   filecont = std::regex_replace(filecont, std::regex("\\$draw_object"), json);

   TString dump_name;
   if ((draw_kind != "draw") && (args.GetBrowserKind() == RWebDisplayArgs::kChrome)) {
      dump_name = "canvasdump";
      FILE *df = gSystem->TempFileName(dump_name);
      if (!df) {
         R__LOG_ERROR(WebGUILog()) << "Fail to create temporary file for dump-dom";
         return false;
      }
      fputs("placeholder", df);
      fclose(df);
   }

   // When true, place HTML file into home directory
   // Some Chrome installation do not allow run html code from files, created in /tmp directory
   static bool chrome_tmp_workaround = false;

   TString tmp_name, html_name;

try_again:

   if ((args.GetBrowserKind() == RWebDisplayArgs::kCEF) || (args.GetBrowserKind() == RWebDisplayArgs::kQt5)) {
      args.SetUrl(""s);
      args.SetPageContent(filecont);

      tmp_name.Clear();
      html_name.Clear();

      R__LOG_DEBUG(0, WebGUILog()) << "Using file content_len " << filecont.length() << " to produce batch image " << fname;

   } else {
      tmp_name = "canvasbody";
      FILE *hf = gSystem->TempFileName(tmp_name);
      if (!hf) {
         R__LOG_ERROR(WebGUILog()) << "Fail to create temporary file for batch job";
         return false;
      }
      fputs(filecont.c_str(), hf);
      fclose(hf);

      html_name = tmp_name + ".html";

      if (chrome_tmp_workaround) {
         std::string homedir = gSystem->GetHomeDirectory();
         auto pos = html_name.Last('/');
         if (pos == kNPOS)
            html_name = TString::Format("/random%d.html", gRandom->Integer(1000000));
         else
            html_name.Remove(0, pos);
         html_name = homedir + html_name.Data();
         gSystem->Unlink(html_name.Data());
         gSystem->Unlink(tmp_name.Data());

         std::ofstream ofs(html_name.Data(), std::ofstream::out);
         ofs << filecont;
      } else {
         if (gSystem->Rename(tmp_name.Data(), html_name.Data()) != 0) {
            R__LOG_ERROR(WebGUILog()) << "Fail to rename temp file " << tmp_name << " into " << html_name;
            gSystem->Unlink(tmp_name.Data());
            return false;
         }
      }

      args.SetUrl("file://"s + gSystem->UnixPathName(html_name.Data()));
      args.SetPageContent(""s);

      R__LOG_DEBUG(0, WebGUILog()) << "Using " << html_name << " content_len " << filecont.length() << " to produce batch image " << fname;
   }

   TString tgtfilename = fname.c_str();
   if (!gSystem->IsAbsoluteFileName(tgtfilename.Data()))
      gSystem->PrependPathName(gSystem->WorkingDirectory(), tgtfilename);

   TString wait_file_name;

   args.SetStandalone(true);
   args.SetHeadless(true);
   args.SetSize(width, height);

   if (draw_kind == "draw") {

      wait_file_name = tgtfilename;

      if (EndsWith(".pdf"))
         args.SetExtraArgs("--print-to-pdf="s + gSystem->UnixPathName(tgtfilename.Data()));
      else
         args.SetExtraArgs("--screenshot="s + gSystem->UnixPathName(tgtfilename.Data()));

   } else if (args.GetBrowserKind() == RWebDisplayArgs::kChrome) {
      // require temporary output file
      args.SetExtraArgs("--dump-dom");
      args.SetRedirectOutput(dump_name.Data());

      // wait_file_name = dump_name;

      gSystem->Unlink(dump_name.Data());
   }

   // remove target image file - we use it as detection when chrome is ready
   gSystem->Unlink(tgtfilename.Data());

   auto handle = RWebDisplayHandle::Display(args);

   if (!handle) {
      R__LOG_DEBUG(0, WebGUILog()) << "Cannot start " << args.GetBrowserName() << " to produce image " << fname;
      return false;
   }

   // delete temporary HTML file
   if (html_name.Length() > 0)
      gSystem->Unlink(html_name.Data());

   if (!wait_file_name.IsNull() && gSystem->AccessPathName(wait_file_name.Data())) {
      R__LOG_ERROR(WebGUILog()) << "Fail to produce image " << fname;
      return false;
   }

   if (draw_kind != "draw") {

      auto dumpcont = handle->GetContent();

      if ((dumpcont.length() > 20) && (dumpcont.length() < 60) && !chrome_tmp_workaround && (args.GetBrowserKind() == RWebDisplayArgs::kChrome)) {
         // chrome creates dummy html file with mostly no content
         // problem running chrome from /tmp directory, lets try work from home directory
         chrome_tmp_workaround = true;
         goto try_again;
      }

      if (dumpcont.length() < 100) {
         R__LOG_ERROR(WebGUILog()) << "Fail to dump HTML code into " << (dump_name.IsNull() ? "CEF" : dump_name.Data());
         return false;
      }

      if (draw_kind == "svg") {
         auto p1 = dumpcont.find("<svg");
         auto p2 = dumpcont.rfind("</svg>");

         std::ofstream ofs(tgtfilename);
         if ((p1 != std::string::npos) && (p2 != std::string::npos) && (p1 < p2)) {
            ofs << dumpcont.substr(p1,p2-p1+6);
         } else {
            R__LOG_ERROR(WebGUILog()) << "Fail to extract SVG from HTML dump " << dump_name;
            ofs << "Failure!!!\n" << dumpcont;
            return false;
         }
      } else {

         auto p1 = dumpcont.rfind(";base64,");
         auto p2 = dumpcont.rfind("></div>");

         if ((p1 != std::string::npos) && (p2 != std::string::npos) && (p1 < p2)) {

            auto base64 = dumpcont.substr(p1+8, p2-p1-9);
            auto binary = TBase64::Decode(base64.c_str());

            std::ofstream ofs(tgtfilename, std::ios::binary);
            ofs.write(binary.Data(), binary.Length());
         } else {
            R__LOG_ERROR(WebGUILog()) << "Fail to extract image from dump HTML code " << dump_name;

            return false;
         }
      }
   }

   R__LOG_DEBUG(0, WebGUILog()) << "Create file " << fname;

   return true;
}

