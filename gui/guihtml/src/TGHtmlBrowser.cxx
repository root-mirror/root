// @(#)root/guitml:$Id$
// Author: Bertrand Bellenot   26/09/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TSystem.h"
#include "TGMenu.h"
#include "TGComboBox.h"
#include "TGButton.h"
#include "TGTextBuffer.h"
#include "TGTextEntry.h"
#include "TGStatusBar.h"
#include "TGFileDialog.h"
#include "TFile.h"
#include "TBrowser.h"
#include "TGHtml.h"
#include "TString.h"
#include "TUrl.h"
#include "TSocket.h"
#include "TGHtmlBrowser.h"
#include "TGText.h"
#include "TError.h"
#include "TVirtualX.h"
#include "snprintf.h"
#ifdef R__SSL
#include "TSSLSocket.h"
#endif
#ifdef WIN32
#include "TWin32SplashThread.h"
#endif

#include <cstdlib>

//_____________________________________________________________________________
//
// TGHtmlBrowser
//
// A very simple HTML browser.
//_____________________________________________________________________________

ClassImp(TGHtmlBrowser);

enum EMyMessageTypes {
   kM_FILE_OPEN,
   kM_FILE_SAVEAS,
   kM_FILE_BROWSE,
   kM_FILE_EXIT,
   kM_FAVORITES_ADD,
   kM_TOOLS_CLEARHIST,
   kM_HELP_ABOUT
};

static const char *gHtmlFTypes[] = {
   "HTML files",    "*.htm*",
   "All files",     "*",
    0,               0
};

const char *HtmlError[] = {
"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3c.org/TR/1999/REC-html401-19991224/loose.dtd\"> ",
"<HTML><HEAD><TITLE>RHTML cannot display the webpage</TITLE> ",
"<META http-equiv=Content-Type content=\"text/html; charset=UTF-8\"></HEAD> ",
"<BODY> ",
"<TABLE cellSpacing=0 cellPadding=0 width=730 border=0> ",
"  <TBODY> ",
"  <TR> ",
"    <TD id=infoIconAlign vAlign=top align=left width=60 rowSpan=2> ",
"    <IMG src=\"info.gif\"> ",
"    </TD> ",
"    <TD id=mainTitleAlign vAlign=center align=left width=*> ",
"      <H1 id=mainTitle>RHTML cannot display the webpage</H1></TD></TR> ",
"  <TR> ",
"    <TD class=errorCodeAndDivider id=errorCodeAlign align=right>&nbsp;  ",
"      <DIV class=divider></DIV></TD></TR> ",
"  <TR> ",
"      <UL> ",
"      </UL> ",
"    <TD>&nbsp; </TD> ",
"    <TD id=MostLikelyAlign vAlign=top align=left> ",
"      <H3 id=likelyCauses>Most likely causes:</H3> ",
"      <UL> ",
"        <LI id=causeNotConnected>You are not connected to the Internet.  ",
"        <LI id=causeSiteProblem>The website is encountering problems.  ",
"        <LI id=causeErrorInAddress>There might be a typing error in the address.  ",
"        <LI id=causeOtherError>  ",
"        </LI></UL></TD></TR> ",
"  <TR> ",
"    <TD id=infoBlockAlign vAlign=top align=right>&nbsp; </TD> ",
"    <TD id=moreInformationAlign vAlign=center align=left> ",
"      <H4> ",
"      <TABLE> ",
"        <TBODY> ",
"        <TR> ",
"          <TD vAlign=top><SPAN id=moreInfoContainer></SPAN><ID  ",
"            id=moreInformation>More information</ID> ",
"      </TD></TR></TBODY></TABLE></H4> ",
"      <DIV class=infoBlock id=infoBlockID> ",
"      <P><ID id=errorExpl1>This problem can be caused by a variety of issues,  ",
"      including:</ID>  ",
"      <UL> ",
"        <LI id=errorExpl2>Internet connectivity has been lost.  ",
"        <LI id=errorExpl3>The website is temporarily unavailable.  ",
"        <LI id=errorExpl4>The Domain Name Server (DNS) is not reachable.  ",
"        <LI id=errorExpl5>The Domain Name Server (DNS) does not have a listing  ",
"        for the website's domain.  ",
"      <P></P> ",
"      <P></P></DIV></TD></TR></TBODY></TABLE></BODY></HTML> ",
0
};

////////////////////////////////////////////////////////////////////////////////
/// TGHtmlBrowser constructor.

TGHtmlBrowser::TGHtmlBrowser(const char *filename, const TGWindow *p, UInt_t w, UInt_t h)
             : TGMainFrame(p, w, h)
{
   SetCleanup(kDeepCleanup);
   fNbFavorites = 1000;
   fMenuBar = new TGMenuBar(this, 35, 50, kHorizontalFrame);

   fMenuFile = new TGPopupMenu(gClient->GetDefaultRoot());
   fMenuFile->AddEntry(" &Open...\tCtrl+O", kM_FILE_OPEN, 0,
                       gClient->GetPicture("ed_open.png"));
   fMenuFile->AddEntry(" Save &As...\tCtrl+A", kM_FILE_SAVEAS, 0,
                       gClient->GetPicture("ed_save.png"));
   fMenuFile->AddEntry(" &Browse...\tCtrl+B", kM_FILE_BROWSE);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry(" E&xit\tCtrl+Q", kM_FILE_EXIT, 0,
                       gClient->GetPicture("bld_exit.png"));
   fMenuFile->Associate(this);

   fMenuFavorites = new TGPopupMenu(gClient->GetDefaultRoot());
   fMenuFavorites->AddEntry("&Add to Favorites", kM_FAVORITES_ADD, 0,
                            gClient->GetPicture("bld_plus.png"));
   fMenuFavorites->AddSeparator();
   fMenuFavorites->AddEntry("http://root.cern.ch", fNbFavorites++, 0,
                            gClient->GetPicture("htmlfile.gif"));
   fMenuFavorites->Associate(this);

   fMenuTools = new TGPopupMenu(gClient->GetDefaultRoot());
   fMenuTools->AddEntry("&Clear History", kM_TOOLS_CLEARHIST, 0,
                        gClient->GetPicture("ed_delete.png"));
   fMenuTools->Associate(this);

   fMenuHelp = new TGPopupMenu(gClient->GetDefaultRoot());
   fMenuHelp->AddEntry(" &About...", kM_HELP_ABOUT, 0, gClient->GetPicture("about.xpm"));
   fMenuHelp->Associate(this);

   fMenuBar->AddPopup("&File", fMenuFile, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Favorites", fMenuFavorites, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Tools", fMenuTools, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Help", fMenuHelp, new TGLayoutHints(kLHintsTop | kLHintsRight));

   AddFrame(fMenuBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   // vertical frame
   fVerticalFrame = new TGVerticalFrame(this,727,600,kVerticalFrame);

   fHorizontalFrame = new TGHorizontalFrame(fVerticalFrame,727,600);

   fBack = new TGPictureButton(fHorizontalFrame,gClient->GetPicture("GoBack.gif"));
   fBack->SetStyle(gClient->GetStyle());
   fBack->SetToolTipText("Go Back");
   fHorizontalFrame->AddFrame(fBack, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));
   fBack->Connect("Clicked()", "TGHtmlBrowser", this, "Back()");

   fForward = new TGPictureButton(fHorizontalFrame,gClient->GetPicture("GoForward.gif"));
   fForward->SetStyle(gClient->GetStyle());
   fForward->SetToolTipText("Go Forward");
   fHorizontalFrame->AddFrame(fForward, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));
   fForward->Connect("Clicked()", "TGHtmlBrowser", this, "Forward()");

   fReload = new TGPictureButton(fHorizontalFrame,gClient->GetPicture("ReloadPage.gif"));
   fReload->SetStyle(gClient->GetStyle());
   fReload->SetToolTipText("Reload Page");
   fHorizontalFrame->AddFrame(fReload, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));
   fReload->Connect("Clicked()", "TGHtmlBrowser", this, "Reload()");

   fStop = new TGPictureButton(fHorizontalFrame,gClient->GetPicture("StopLoading.gif"));
   fStop->SetStyle(gClient->GetStyle());
   fStop->SetToolTipText("Stop Loading");
   fHorizontalFrame->AddFrame(fStop, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));
   fStop->Connect("Clicked()", "TGHtmlBrowser", this, "Stop()");

   fHome = new TGPictureButton(fHorizontalFrame,gClient->GetPicture("GoHome.gif"));
   fHome->SetStyle(gClient->GetStyle());
   fHome->SetToolTipText("Go to ROOT HomePage\n  (http://root.cern.ch)");
   fHorizontalFrame->AddFrame(fHome, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));
   fHome->Connect("Clicked()", "TGHtmlBrowser", this, "Selected(=\"http://root.cern.ch\")");

   // combo box
   fURLBuf   = new TGTextBuffer(256);
   fComboBox = new TGComboBox(fHorizontalFrame, "");
   fURL      = fComboBox->GetTextEntry();
   fURLBuf   = fURL->GetBuffer();
   fComboBox->Resize(200, fURL->GetDefaultHeight());
   fURL->Connect("ReturnPressed()", "TGHtmlBrowser", this, "URLChanged()");

   if (filename) {
      fComboBox->AddEntry(filename, 1);
      fURL->SetText(filename);
   }
   fComboBox->Select(0);
   fComboBox->Connect("Selected(char *)", "TGHtmlBrowser", this, "Selected(char *)");

   fHorizontalFrame->AddFrame(fComboBox, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

   fVerticalFrame->AddFrame(fHorizontalFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));

   // embedded canvas
   fHtml = new TGHtml(fVerticalFrame, 10, 10, -1);
   fVerticalFrame->AddFrame(fHtml, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));

   AddFrame(fVerticalFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));

   // status bar
   fStatusBar = new TGStatusBar(this,100,20);
   Int_t partsusBar[] = {75,25};
   fStatusBar->SetParts(partsusBar,2);
   AddFrame(fStatusBar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX));

   fHtml->Connect("MouseOver(const char *)", "TGHtmlBrowser", this, "MouseOver(const char *)");
   fHtml->Connect("MouseDown(const char *)", "TGHtmlBrowser", this, "MouseDown(const char *)");

   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();
   Resize(w, h);

   if (filename)
      Selected(filename);
}

////////////////////////////////////////////////////////////////////////////////
/// Read (open) remote files.

Ssiz_t ReadSize(const char *url)
{
   char buf[4096];
   TUrl fUrl(url);

   // Give full URL so Apache's virtual hosts solution works.
   TString msg = "HEAD ";
   msg += fUrl.GetProtocol();
   msg += "://";
   msg += fUrl.GetHost();
   msg += ":";
   msg += fUrl.GetPort();
   msg += "/";
   msg += fUrl.GetFile();
   msg += " HTTP/1.0";
   msg += "\r\n";
   msg += "User-Agent: ROOT-TWebFile/1.1";
   msg += "\r\n\r\n";

   TSocket *s;
   TString uri(url);
   if (!uri.BeginsWith("http://") && !uri.BeginsWith("https://"))
      return 0;
   if (uri.BeginsWith("https://")) {
#ifdef R__SSL
      s = new TSSLSocket(fUrl.GetHost(), fUrl.GetPort());
#else
      ::Error("ReadSize", "library compiled without SSL, https not supported");
      return 0;
#endif
   }
   else {
      s = new TSocket(fUrl.GetHost(), fUrl.GetPort());
   }
   if (!s->IsValid()) {
      delete s;
      return 0;
   }
   if (s->SendRaw(msg.Data(), msg.Length()) == -1) {
      delete s;
      return 0;
   }
   if (s->RecvRaw(buf, 4096) == -1) {
      delete s;
      return 0;
   }
   TString reply(buf);
   Ssiz_t idx = reply.Index("Content-length:", 0, TString::kIgnoreCase);
   if (idx > 0) {
      idx += 15;
      TString slen = reply(idx, reply.Length() - idx);
      delete s;
      return (Ssiz_t)atol(slen.Data());
   }
   delete s;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read (open) remote files.

static char *ReadRemote(const char *url)
{
   static char *buf = 0;
   TUrl fUrl(url);

   Ssiz_t size = ReadSize(url);
   if (size <= 0) size = 1024*1024;

   TString msg = "GET ";
   msg += fUrl.GetProtocol();
   msg += "://";
   msg += fUrl.GetHost();
   msg += ":";
   msg += fUrl.GetPort();
   msg += "/";
   msg += fUrl.GetFile();
   msg += "\r\n";

   TSocket *s;
   TString uri(url);
   if (!uri.BeginsWith("http://") && !uri.BeginsWith("https://"))
      return 0;
   if (uri.BeginsWith("https://")) {
#ifdef R__SSL
      s = new TSSLSocket(fUrl.GetHost(), fUrl.GetPort());
#else
      ::Error("ReadRemote", "library compiled without SSL, https not supported");
     return 0;
#endif
   }
   else {
      s = new TSocket(fUrl.GetHost(), fUrl.GetPort());
   }
   if (!s->IsValid()) {
      delete s;
      return 0;
   }
   if (s->SendRaw(msg.Data(), msg.Length()) == -1) {
      delete s;
      return 0;
   }
   buf = (char *)calloc(size+1, sizeof(char));
   if (s->RecvRaw(buf, size) == -1) {
      free(buf);
      delete s;
      return 0;
   }
   delete s;
   return buf;
}

////////////////////////////////////////////////////////////////////////////////
/// Open (browse) selected URL.

void TGHtmlBrowser::Selected(const char *uri)
{
   char *buf = 0;
   FILE *f;

   if (CheckAnchors(uri))
      return;

   TString surl(gSystem->UnixPathName(uri));
   if (!surl.BeginsWith("http://") && !surl.BeginsWith("https://") &&
       !surl.BeginsWith("ftp://") && !surl.BeginsWith("file://")) {
      if (surl.BeginsWith("file:"))
         surl.ReplaceAll("file:", "file://");
      else
         surl.Prepend("file://");
   }
   if (surl.EndsWith(".root")) {
      // in case of root file, just open it and refresh browsers
      gVirtualX->SetCursor(fHtml->GetId(), gVirtualX->CreateCursor(kWatch));
      gROOT->ProcessLine(Form("TFile::Open(\"%s\");", surl.Data()));
      Clicked((char *)surl.Data());
      gROOT->RefreshBrowsers();
      gVirtualX->SetCursor(fHtml->GetId(), gVirtualX->CreateCursor(kPointer));
      return;
   }
   gVirtualX->SetCursor(fHtml->GetId(), gVirtualX->CreateCursor(kWatch));
   TUrl url(surl.Data());
   if (surl.EndsWith(".pdf", TString::kIgnoreCase)) {
      // special case: open pdf files with external viewer
      // works only on Windows for the time being...
      if (!gVirtualX->InheritsFrom("TGX11")) {
         TString cmd = TString::Format("explorer %s", surl.Data());
         gSystem->Exec(cmd.Data());
      }
      gVirtualX->SetCursor(fHtml->GetId(), gVirtualX->CreateCursor(kPointer));
      return;
   }
   if (surl.EndsWith(".gif") || surl.EndsWith(".jpg") || surl.EndsWith(".png")) {
      // special case: single picture
      fHtml->Clear();
      char imgHtml[1024];
      snprintf(imgHtml, 1000, "<IMG src=\"%s\"> ", surl.Data());
      fHtml->ParseText(imgHtml);
      fHtml->SetBaseUri(url.GetUrl());
      fURL->SetText(surl.Data());
      if (!fComboBox->FindEntry(surl.Data()))
         fComboBox->AddEntry(surl.Data(), fComboBox->GetNumberOfEntries()+1);
      fHtml->Layout();
      gVirtualX->SetCursor(fHtml->GetId(), gVirtualX->CreateCursor(kPointer));
      return;
   }
   if (!strcmp(url.GetProtocol(), "http") ||
       !strcmp(url.GetProtocol(), "https")) {
      // standard web page
      buf = ReadRemote(url.GetUrl());
      if (buf) {
         fHtml->Clear();
         fHtml->Layout();
         fHtml->SetBaseUri(url.GetUrl());
         fHtml->ParseText(buf);
         free(buf);
         fURL->SetText(surl.Data());
         if (!fComboBox->FindEntry(surl.Data()))
            fComboBox->AddEntry(surl.Data(), fComboBox->GetNumberOfEntries()+1);
      }
      else {
         fHtml->Clear();
         fHtml->Layout();
         fHtml->SetBaseUri("");
         for (int i=0; HtmlError[i]; i++) {
            fHtml->ParseText((char *)HtmlError[i]);
         }
      }
   }
   else {
      // local file
      f = fopen(url.GetFile(), "r");
      if (f) {
         TString fpath = url.GetUrl();
         fpath.ReplaceAll(gSystem->BaseName(fpath.Data()), "");
         fpath.ReplaceAll("file://", "");
         fHtml->Clear();
         fHtml->Layout();
         fHtml->SetBaseUri(fpath.Data());
         buf = (char *)calloc(4096, sizeof(char));
         if (buf) {
            while (fgets(buf, 4096, f)) {
               fHtml->ParseText(buf);
            }
            free(buf);
         }
         fclose(f);
         fURL->SetText(surl.Data());
         if (!fComboBox->FindEntry(surl.Data()))
            fComboBox->AddEntry(surl.Data(), fComboBox->GetNumberOfEntries()+1);
      }
      else {
         fHtml->Clear();
         fHtml->Layout();
         fHtml->SetBaseUri("");
         for (int i=0; HtmlError[i]; i++) {
            fHtml->ParseText((char *)HtmlError[i]);
         }
      }
   }
   // restore cursor
   gVirtualX->SetCursor(fHtml->GetId(), gVirtualX->CreateCursor(kPointer));
   fHtml->Layout();
   Ssiz_t idx = surl.Last('#');
   if (idx > 0) {
      idx +=1; // skip #
      TString anchor = surl(idx, surl.Length() - idx);
      fHtml->GotoAnchor(anchor.Data());
   }
   SetWindowName(Form("%s - RHTML",surl.Data()));
}

////////////////////////////////////////////////////////////////////////////////
/// URL combobox has changed.

void TGHtmlBrowser::URLChanged()
{
   const char *string = fURL->GetText();
   if (string) {
      Selected(gSystem->UnixPathName(string));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle "Back" navigation button.

void TGHtmlBrowser::Back()
{
   Int_t index = 0;
   const char *string = fURL->GetText();
   TGLBEntry * lbe1 = fComboBox->FindEntry(string);
   if (lbe1)
      index = lbe1->EntryId();
   if (index > 0) {
      fComboBox->Select(index - 1, kTRUE);
      TGTextLBEntry *entry = (TGTextLBEntry *)fComboBox->GetSelectedEntry();
      if (entry) {
         string = entry->GetTitle();
         if (string)
            Selected(string);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if we just change position in the page (using anchor)
/// and return kTRUE if any anchor has been found and followed.

Bool_t TGHtmlBrowser::CheckAnchors(const char *uri)
{
   TString surl(gSystem->UnixPathName(uri));

   if (!fHtml->GetBaseUri())
      return kFALSE;
   TString actual = fHtml->GetBaseUri();
   Ssiz_t idx = surl.Last('#');
   Ssiz_t idy = actual.Last('#');
   TString short1(surl.Data());
   TString short2(actual.Data());
   if (idx > 0)
      short1 = surl(0, idx);
   if (idy > 0)
      short2 = actual(0, idy);

   if (short1 == short2) {
      if (idx > 0) {
         idx +=1; // skip #
         TString anchor = surl(idx, surl.Length() - idx);
         fHtml->GotoAnchor(anchor.Data());
      }
      else {
         fHtml->ScrollToPosition(TGLongPosition(0, 0));
      }
      fHtml->SetBaseUri(surl.Data());
      if (!fComboBox->FindEntry(surl.Data()))
         fComboBox->AddEntry(surl.Data(), fComboBox->GetNumberOfEntries()+1);
      fURL->SetText(surl.Data());
      fComboBox->Select(fComboBox->GetNumberOfEntries(), kFALSE);
      SetWindowName(Form("%s - RHTML",surl.Data()));
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle "Forward" navigation button.

void TGHtmlBrowser::Forward()
{
   Int_t index = 0;
   const char *string = fURL->GetText();
   TGLBEntry * lbe1 = fComboBox->FindEntry(string);
   if (lbe1)
      index = lbe1->EntryId();
   if (index < fComboBox->GetNumberOfEntries()) {
      fComboBox->Select(index + 1, kTRUE);
      TGTextLBEntry *entry = (TGTextLBEntry *)fComboBox->GetSelectedEntry();
      if (entry) {
         string = entry->GetTitle();
         if (string)
            Selected(string);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle "Reload" navigation button.

void TGHtmlBrowser::Reload()
{
   const char *string = fURL->GetText();
   if (string)
      Selected(string);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle "Reload" navigation button.

void TGHtmlBrowser::Stop()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Handle "MouseOver" TGHtml signal.

void TGHtmlBrowser::MouseOver(const char *url)
{
   fStatusBar->SetText(url, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle "MouseDown" TGHtml signal.

void TGHtmlBrowser::MouseDown(const char *url)
{
   Selected(url);
}

////////////////////////////////////////////////////////////////////////////////
/// Process Events.

Bool_t TGHtmlBrowser::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   switch (GET_MSG(msg)) {
   case kC_COMMAND:
      {
         switch (GET_SUBMSG(msg)) {

            case kCM_MENU:
            case kCM_BUTTON:

               switch(parm1) {

                  case kM_FILE_EXIT:
                     CloseWindow();
                     break;

                  case kM_FILE_OPEN:
                     {
                        static TString dir(".");
                        TGFileInfo fi;
                        fi.fFileTypes = gHtmlFTypes;
                        fi.SetIniDir(dir);
                        new TGFileDialog(fClient->GetRoot(), this,
                                         kFDOpen, &fi);
                        dir = fi.fIniDir;
                        if (fi.fFilename) {
                           Selected(Form("file://%s",
                              gSystem->UnixPathName(fi.fFilename)));
                        }
                     }
                     break;

                  case kM_FILE_SAVEAS:
                     {
                        static TString sdir(".");
                        TGFileInfo fi;
                        fi.fFileTypes = gHtmlFTypes;
                        fi.SetIniDir(sdir);
                        new TGFileDialog(fClient->GetRoot(), this,
                                         kFDSave, &fi);
                        sdir = fi.fIniDir;
                        if (fi.fFilename) {
                           TGText txt(fHtml->GetText());
                           txt.Save(gSystem->UnixPathName(fi.fFilename));
                        }
                     }
                     break;

                  case kM_FAVORITES_ADD:
                     fMenuFavorites->AddEntry(Form("%s",
                           fURL->GetText()), fNbFavorites++, 0,
                           gClient->GetPicture("htmlfile.gif"));
                     break;

                  case kM_TOOLS_CLEARHIST:
                     fComboBox->RemoveEntries(1,fComboBox->GetNumberOfEntries());
                     break;

                  case kM_FILE_BROWSE:
                     new TBrowser();
                     break;

                  case kM_HELP_ABOUT:
                     {
#ifdef R__UNIX
                        TString rootx = TROOT::GetBinDir() + "/root -a &";
                        gSystem->Exec(rootx);
#else
#ifdef WIN32
                        new TWin32SplashThread(kTRUE);
#else
                        char str[32];
                        snprintf(str, 32, "About ROOT %s...", gROOT->GetVersion());
                        TRootHelpDialog *hd = new TRootHelpDialog(this, str,
                                                                  600, 400);
                        hd->SetText(gHelpAbout);
                        hd->Popup();
#endif
#endif
                     }
                     break;

                  default:
                     {
                        if (parm1 < 1000) break;
                        TGMenuEntry *entry = fMenuFavorites->GetEntry(parm1);
                        if (!entry) break;
                        const char *shortcut = entry->GetName();
                        if (shortcut)
                           Selected(shortcut);
                     }
                     break;
               }
               break;
         }
         break;
      }
   }
   return kTRUE;
}

