// @(#)root/gui:$Id$
// Author: Bertrand Bellenot   26/09/2007

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootBrowser                                                         //
//                                                                      //
// This class creates a ROOT object browser, constitued by three main   //
// tabs.                                                                //
//                                                                      //
// All tabs can 'swallow' frames, thanks to the new method:             //
//   ExecPlugin(const char *name = 0, const char *fname = 0,            //
//              const char *cmd = 0, Int_t pos = kRight,                //
//              Int_t subpos = -1)                                      //
// allowing to select plugins (can be a macro or a command)             //
// to be executed, and where to embed the frame created by              //
// the plugin (tab and tab element). Examples:                          //
//                                                                      //
// create a new browser:                                                //
// TBrowser b;                                                          //
//                                                                      //
// create a new TCanvas in a new top right tab element:                 //
// b.ExecPlugin("Canvas", 0, "new TCanvas()");                          //
//                                                                      //
// create a new top right tab element embedding the                     //
// TGMainFrame created by the macro 'myMacro.C':                        //
// b.ExecPlugin("MyPlugin", "myMacro.C");                               //
//                                                                      //
// create a new bottom tab element embedding the                        //
// TGMainFrame created by the macro 'myMacro.C':                        //
// b.ExecPlugin("MyPlugin", "myMacro.C", 0, TRootBrowser::kBottom);     //
//                                                                      //
// this browser implementation can be selected via the env              //
// 'Browser.Name' in .rootrc, (TRootBrowser or TRootBrowserLite)        //
// the default being TRootBrowserLite (old browser)                     //
// a list of options (plugins) for the new TRootBrowser is also         //
// specified via the env 'Browser.Options' in .rootrc, the default      //
// being: FECI                                                          //
// Here is the list of available options:                               //
// F: File browser E: Text Editor H: HTML browser C: Canvas I: I/O      //
// redirection P: Proof G: GL viewer                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TBrowser.h"
#include "TClass.h"
#include "TGClient.h"
#include "TGFrame.h"
#include "TGTab.h"
#include "TGMenu.h"
#include "TGLayout.h"
#include "TGSplitter.h"
#include "TGStatusBar.h"
#include "Varargs.h"
#include "TInterpreter.h"
#include "TGFileDialog.h"
#include "TObjString.h"
#include "TVirtualPad.h"
#include "TEnv.h"
#include <KeySymbols.h>

#include "RConfigure.h"

#include "TRootBrowser.h"
#include "TGFileBrowser.h"
#include "TGInputDialog.h"
#include "TRootHelpDialog.h"
#include "TVirtualPadEditor.h"
#include "HelpText.h"
#include "Getline.h"
#include "TVirtualX.h"
#include "strlcpy.h"
#include "snprintf.h"

#ifdef WIN32
#include <TWin32SplashThread.h>
#endif

static const char *gOpenFileTypes[] = {
   "ROOT files",   "*.root",
   "All files",    "*",
   0,              0
};

static const char *gPluginFileTypes[] = {
   "ROOT files",   "*.C",
   "All files",    "*",
   0,              0
};

//_____________________________________________________________________________
//
// TRootBrowser
//
// The main ROOT object browser.
//_____________________________________________________________________________

ClassImp(TRootBrowser);

////////////////////////////////////////////////////////////////////////////////
/// Create browser with a specified width and height.

TRootBrowser::TRootBrowser(TBrowser *b, const char *name, UInt_t width,
                           UInt_t height, Option_t *opt, Bool_t initshow) :
   TGMainFrame(gClient->GetDefaultRoot(), width, height), TBrowserImp(b)
{
   fShowCloseTab = kTRUE;
   fActBrowser = 0;
   fIconPic = 0;
   CreateBrowser(name);
   Resize(width, height);
   if (initshow) {
      InitPlugins(opt);
      MapWindow();
   }
   TQObject::Connect("TCanvas", "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",
                     "TRootBrowser", this,
                     "EventInfo(Int_t, Int_t, Int_t, TObject*)");
   gVirtualX->SetInputFocus(GetId());
}

////////////////////////////////////////////////////////////////////////////////
/// Create browser with a specified width and height and at position x, y.

TRootBrowser::TRootBrowser(TBrowser *b, const char *name, Int_t x, Int_t y,
                           UInt_t width, UInt_t height, Option_t *opt,
                           Bool_t initshow) :
   TGMainFrame(gClient->GetDefaultRoot(), width, height), TBrowserImp(b)
{
   fShowCloseTab = kTRUE;
   fActBrowser = 0;
   fIconPic = 0;
   CreateBrowser(name);
   MoveResize(x, y, width, height);
   SetWMPosition(x, y);
   if (initshow) {
      InitPlugins(opt);
      MapWindow();
   }
   TQObject::Connect("TCanvas", "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",
                     "TRootBrowser", this,
                     "EventInfo(Int_t, Int_t, Int_t, TObject*)");
   gVirtualX->SetInputFocus(GetId());
}


////////////////////////////////////////////////////////////////////////////////

void TRootBrowser::CreateBrowser(const char *name)
{
   // Create the actual interface.

   fVf = new TGVerticalFrame(this, 100, 100);

   fLH0 = new TGLayoutHints(kLHintsNormal);
   fLH1 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
   fLH2 = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 1, 1);
   fLH3 = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX);
   fLH4 = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2);
   fLH5 = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY);
   fLH6 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX);
   fLH7 = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY);

   // Menubar Frame
   fTopMenuFrame = new TGHorizontalFrame(fVf, 100, 20);

   fPreMenuFrame = new TGHorizontalFrame(fTopMenuFrame, 0, 20, kRaisedFrame);
   fMenuBar   = new TGMenuBar(fPreMenuFrame, 10, 10, kHorizontalFrame);
   fMenuFile  = new TGPopupMenu(gClient->GetDefaultRoot());
   fMenuFile->AddEntry("&Browse...\tCtrl+B", kBrowse);
   fMenuFile->AddEntry("&Open...\tCtrl+O", kOpenFile);
   fMenuFile->AddSeparator();

   fMenuHelp = new TGPopupMenu(fClient->GetRoot());
   fMenuHelp->AddEntry("&About ROOT...",        kHelpAbout);
   fMenuHelp->AddSeparator();
   fMenuHelp->AddEntry("Help On Browser...",    kHelpOnBrowser);
   fMenuHelp->AddEntry("Help On Canvas...",     kHelpOnCanvas);
   fMenuHelp->AddEntry("Help On Menus...",      kHelpOnMenus);
   fMenuHelp->AddEntry("Help On Graphics Editor...", kHelpOnGraphicsEd);
   fMenuHelp->AddEntry("Help On Objects...",    kHelpOnObjects);
   fMenuHelp->AddEntry("Help On PostScript...", kHelpOnPS);
   fMenuHelp->AddEntry("Help On Remote Session...", kHelpOnRemote);
   fMenuFile->AddPopup("Browser Help...", fMenuHelp);

   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("&Clone\tCtrl+N", kClone);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("New &Editor\tCtrl+E", kNewEditor);
   fMenuFile->AddEntry("New &Canvas\tCtrl+C", kNewCanvas);
   fMenuFile->AddEntry("New &HTML\tCtrl+H", kNewHtml);
   fMenuFile->AddSeparator();
   fMenuExecPlugin = new TGPopupMenu(fClient->GetRoot());
   fMenuExecPlugin->AddEntry("&Macro...", kExecPluginMacro);
   fMenuExecPlugin->AddEntry("&Command...", kExecPluginCmd);
   fMenuFile->AddPopup("Execute &Plugin...", fMenuExecPlugin);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("Close &Tab\tCtrl+T", kCloseTab);
   fMenuFile->AddEntry("Close &Window\tCtrl+W", kCloseWindow);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("&Quit Root\tCtrl+Q", kQuitRoot);
   fMenuBar->AddPopup("&Browser", fMenuFile, fLH1);
   fMenuFile->Connect("Activated(Int_t)", "TRootBrowser", this,
                      "HandleMenu(Int_t)");
   fPreMenuFrame->AddFrame(fMenuBar, fLH2);
   fTopMenuFrame->AddFrame(fPreMenuFrame, fLH0);

   if (!TClass::GetClass("TGHtmlBrowser"))
      fMenuFile->DisableEntry(kNewHtml);

   fMenuFrame = new TGHorizontalFrame(fTopMenuFrame, 100, 20, kRaisedFrame);
   fTopMenuFrame->AddFrame(fMenuFrame, fLH5);

   fVf->AddFrame(fTopMenuFrame, fLH3);
   fActMenuBar = fMenuBar;

   // Toolbar Frame
   fToolbarFrame = new TGHorizontalFrame(fVf, 100, 20, kHorizontalFrame |
                                         kRaisedFrame);
   fVf->AddFrame(fToolbarFrame, fLH3);

   fHf = new TGHorizontalFrame(fVf, 100, 100);
   // Tabs & co...
#if defined(R__HAS_COCOA)
   fV1 = new TGVerticalFrame(fHf, 252, 100, kFixedWidth);
#else
   fV1 = new TGVerticalFrame(fHf, 250, 100, kFixedWidth);
#endif
   fV2 = new TGVerticalFrame(fHf, 600, 100);
   fH1 = new TGHorizontalFrame(fV2, 100, 100);
   fH2 = new TGHorizontalFrame(fV2, 100, 100, kFixedHeight);

   // Left tab
   fTabLeft = new TGTab(fV1,100,100);
   //fTabLeft->AddTab("Tab 1");
   fTabLeft->Resize(fTabLeft->GetDefaultSize());
   fV1->AddFrame(fTabLeft, fLH4);

   // Vertical splitter
   fVSplitter = new TGVSplitter(fHf, 4, 4);
   fVSplitter->SetFrame(fV1, kTRUE);
   fHf->AddFrame(fV1, fLH7);
   fHf->AddFrame(fVSplitter, fLH7);

   // Right tab
   fTabRight = new TGTab(fH1, 500, 100);
   //fTabRight->AddTab("Tab 1");
   fTabRight->Resize(fTabRight->GetDefaultSize());
   fH1->AddFrame(fTabRight, fLH5);
   fTabRight->Connect("Selected(Int_t)", "TRootBrowser", this, "DoTab(Int_t)");
   fTabRight->Connect("CloseTab(Int_t)", "TRootBrowser", this, "CloseTab(Int_t)");
   fV2->AddFrame(fH1, fLH4);

   // Horizontal splitter
   fHSplitter = new TGHSplitter(fV2, 4, 4);
   fV2->AddFrame(fHSplitter, fLH3);

   // Bottom tab
   fTabBottom = new TGTab(fH2, 100, 100);
   //fTabBottom->AddTab("Tab 1");
   fH2->AddFrame(fTabBottom, fLH4);
   fV2->AddFrame(fH2, fLH3);

   fHSplitter->SetFrame(fH2, kFALSE);
   fHf->AddFrame(fV2, fLH5);
   fVf->AddFrame(fHf, fLH5);
   AddFrame(fVf, fLH5);

   // status bar
   fStatusBar = new TGStatusBar(this, 400, 20);
   int parts[] = { 33, 10, 10, 47 };
   fStatusBar->SetParts(parts, 4);
   AddFrame(fStatusBar, fLH6);

   fNbInitPlugins = 0;
   fEditFrame = 0;
   fEditTab   = 0;
   fEditPos   = -1;
   fEditSubPos= -1;
   fNbTab[0]  = fNbTab[1] = fNbTab[2] = 0;
   fCrTab[0]  = fCrTab[1] = fCrTab[2] = -1;

   // Set a name to the main frame
   SetWindowName(name);
   SetIconName(name);
   fIconPic = SetIconPixmap("rootdb_s.xpm");
   SetClassHints("ROOT", "Browser");

   if (!strcmp(gROOT->GetDefCanvasName(), "c1"))
      gROOT->SetDefCanvasName("Canvas_1");

   SetWMSizeHints(600, 350, 10000, 10000, 2, 2);
   MapSubwindows();
   Resize(GetDefaultSize());
   AddInput(kKeyPressMask | kKeyReleaseMask);

   fVf->HideFrame(fToolbarFrame);
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up all widgets, frames and layouthints that were used

TRootBrowser::~TRootBrowser()
{
   if (fIconPic) gClient->FreePicture(fIconPic);
   delete fLH0;
   delete fLH1;
   delete fLH2;
   delete fLH3;
   delete fLH4;
   delete fLH5;
   delete fLH6;
   delete fLH7;
   delete fMenuHelp;
   delete fMenuExecPlugin;
   delete fMenuFile;
   delete fMenuBar;
   delete fMenuFrame;
   delete fPreMenuFrame;
   delete fTopMenuFrame;
   delete fToolbarFrame;
   delete fVSplitter;
   delete fHSplitter;
   delete fTabLeft;
   delete fTabRight;
   delete fTabBottom;
   delete fH1;
   delete fH2;
   delete fV1;
   delete fV2;
   delete fHf;
   delete fStatusBar;
   delete fVf;
}

////////////////////////////////////////////////////////////////////////////////
/// Add items to the actual browser. This function has to be called
/// by the Browse() member function of objects when they are
/// called by a browser. If check < 0 (default) no check box is drawn,
/// if 0 then unchecked checkbox is added, if 1 checked checkbox is added.

void TRootBrowser::Add(TObject *obj, const char *name, Int_t check)
{
   if (obj->InheritsFrom("TObjectSpy"))
      return;
   if (fActBrowser)
      fActBrowser->Add(obj, name, check);
}

////////////////////////////////////////////////////////////////////////////////
/// Browse object. This, in turn, will trigger the calling of
/// TRootBrowser::Add() which will fill the IconBox and the tree.
/// Emits signal "BrowseObj(TObject*)".

void TRootBrowser::BrowseObj(TObject *obj)
{
   if (fActBrowser)
      fActBrowser->BrowseObj(obj);
   Emit("BrowseObj(TObject*)", (Long_t)obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Clone the browser. A new Browser will be created, with the same
/// plugins executed in the current one.

void TRootBrowser::CloneBrowser()
{
   Int_t loop = 1;
   TBrowserPlugin *plugin = 0;
   TBrowser *b = new TBrowser();
   TIter next(&fPlugins);
   while ((plugin = (TBrowserPlugin *)next())) {
      if (loop > fNbInitPlugins)
         b->ExecPlugin(plugin->GetName(), "", plugin->fCommand.Data(), plugin->fTab,
                       plugin->fSubTab);
      ++loop;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove tab element id from right tab.

void TRootBrowser::CloseTab(Int_t id)
{
   RemoveTab(kRight, id);
}

////////////////////////////////////////////////////////////////////////////////
/// Properly close the mainframes embedded in the different tabs

void TRootBrowser::CloseTabs()
{
   TGFrameElement *el;
   TGCompositeFrame *container;
   Int_t i;
   TQObject::Disconnect("TCanvas", "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",
                        this, "EventInfo(Int_t, Int_t, Int_t, TObject*)");
   Disconnect(fMenuFile, "Activated(Int_t)", this, "HandleMenu(Int_t)");
   Disconnect(fTabRight, "Selected(Int_t)", this, "DoTab(Int_t)");
   if (fPlugins.IsEmpty()) return;
   fActBrowser = 0;
   for (i=0;i<fTabLeft->GetNumberOfTabs();i++) {
      container = fTabLeft->GetTabContainer(i);
      if (!container) continue;
      el = (TGFrameElement *)container->GetList()->First();
      if (el && el->fFrame) {
         el->fFrame->SetFrameElement(0);
         if (el->fFrame->InheritsFrom("TVirtualPadEditor")) {
            TVirtualPadEditor::Terminate();
         }
         else if (el->fFrame->InheritsFrom("TGMainFrame")) {
            el->fFrame->UnmapWindow();
            ((TGMainFrame *)el->fFrame)->CloseWindow();
            gSystem->ProcessEvents();
         }
         else
            delete el->fFrame;
         el->fFrame = 0;
         if (el->fLayout && (el->fLayout != fgDefaultHints) &&
            (el->fLayout->References() > 0)) {
            el->fLayout->RemoveReference();
            if (!el->fLayout->References()) {
               delete el->fLayout;
            }
         }
         container->GetList()->Remove(el);
         delete el;
      }
   }
   for (i=0;i<fTabRight->GetNumberOfTabs();i++) {
      container = fTabRight->GetTabContainer(i);
      if (!container) continue;
      el = (TGFrameElement *)container->GetList()->First();
      if (el && el->fFrame) {
         el->fFrame->SetFrameElement(0);
         if (el->fFrame->InheritsFrom("TGMainFrame")) {
            el->fFrame->UnmapWindow();
            Bool_t sleep = (el->fFrame->InheritsFrom("TRootCanvas")) ? kTRUE : kFALSE;
            if (sleep)
               gSystem->Sleep(150);
            ((TGMainFrame *)el->fFrame)->CloseWindow();
            if (sleep)
               gSystem->Sleep(150);
            gSystem->ProcessEvents();
         }
         else
            delete el->fFrame;
         el->fFrame = 0;
         if (el->fLayout && (el->fLayout != fgDefaultHints) &&
            (el->fLayout->References() > 0)) {
            el->fLayout->RemoveReference();
            if (!el->fLayout->References()) {
               delete el->fLayout;
            }
         }
         container->GetList()->Remove(el);
         delete el;
      }
   }
   for (i=0;i<fTabBottom->GetNumberOfTabs();i++) {
      container = fTabBottom->GetTabContainer(i);
      if (!container) continue;
      el = (TGFrameElement *)container->GetList()->First();
      if (el && el->fFrame) {
         el->fFrame->SetFrameElement(0);
         if (el->fFrame->InheritsFrom("TGMainFrame")) {
            el->fFrame->UnmapWindow();
            ((TGMainFrame *)el->fFrame)->CloseWindow();
            gSystem->Sleep(150);
            gSystem->ProcessEvents();
         }
         else
            delete el->fFrame;
         el->fFrame = 0;
         if (el->fLayout && (el->fLayout != fgDefaultHints) &&
            (el->fLayout->References() > 0)) {
            el->fLayout->RemoveReference();
            if (!el->fLayout->References()) {
               delete el->fLayout;
            }
         }
         container->GetList()->Remove(el);
         delete el;
      }
   }
   fPlugins.Delete();
   Emit("CloseWindow()");
}

////////////////////////////////////////////////////////////////////////////////
/// Called when window is closed via the window manager.

void TRootBrowser::CloseWindow()
{
   TQObject::Disconnect("TCanvas", "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",
                        this, "EventInfo(Int_t, Int_t, Int_t, TObject*)");
   CloseTabs();
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle Tab navigation.

void TRootBrowser::DoTab(Int_t id)
{
   TGTab *sender = (TGTab *)gTQSender;
   if ((sender) && (sender == fTabRight)) {
      SwitchMenus(sender->GetTabContainer(id));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Display a tooltip with infos about the primitive below the cursor.

void TRootBrowser::EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected)
{
   const Int_t kTMAX=256;
   static char atext[kTMAX];
   if (selected == 0 || event == kMouseLeave) {
      SetStatusText("", 0);
      SetStatusText("", 1);
      SetStatusText("", 2);
      SetStatusText("", 3);
      return;
   }
   SetStatusText(selected->GetTitle(), 0);
   SetStatusText(selected->GetName(), 1);
   if (event == kKeyPress)
      snprintf(atext, kTMAX, "%c", (char) px);
   else
      snprintf(atext, kTMAX, "%d,%d", px, py);
   SetStatusText(atext, 2);
   SetStatusText(selected->GetObjectInfo(px,py), 3);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute a macro and embed the created frame in the tab "pos"
/// and tab element "subpos".

Long_t TRootBrowser::ExecPlugin(const char *name, const char *fname,
                                const char *cmd, Int_t pos, Int_t subpos)
{
   Long_t retval = 0;
   TBrowserPlugin *p;
   TString command, pname;
   if (cmd && strlen(cmd)) {
      command = cmd;
      if (name) pname = name;
      else pname = TString::Format("Plugin %d", fPlugins.GetSize());
      p = new TBrowserPlugin(pname.Data(), command.Data(), pos, subpos);
   }
   else if (fname && strlen(fname)) {
      pname = name ? name : gSystem->BaseName(fname);
      Ssiz_t t = pname.Last('.');
      if (t > 0) pname.Remove(t);
      command.Form("gROOT->Macro(\"%s\");", gSystem->UnixPathName(fname));
      p = new TBrowserPlugin(pname.Data(), command.Data(), pos, subpos);
   }
   else return 0;
   if (IsWebGUI() && command.Contains("new TCanvas"))
      return gROOT->ProcessLine(command.Data());
   StartEmbedding(pos, subpos);
   fPlugins.Add(p);
   retval = gROOT->ProcessLine(command.Data());
   if (command.Contains("new TCanvas")) {
      pname = gPad->GetName();
      p->SetName(pname.Data());
   }
   SetTabTitle(pname.Data(), pos, subpos);
   StopEmbedding();
   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns drawing option.

Option_t *TRootBrowser::GetDrawOption() const
{
   if (fActBrowser)
      return fActBrowser->GetDrawOption();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the TGTab at position pos.

TGTab* TRootBrowser::GetTab(Int_t pos) const
{
   switch (pos) {
      case kLeft:   return fTabLeft;
      case kRight:  return fTabRight;
      case kBottom: return fTabBottom;
      default:      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle keyboard events.

Bool_t TRootBrowser::HandleKey(Event_t *event)
{
   char   input[10];
   UInt_t keysym;

   if (event->fType == kGKeyPress) {
      gVirtualX->LookupString(event, input, sizeof(input), keysym);

      if (!event->fState && (EKeySym)keysym == kKey_F5) {
         Refresh(kTRUE);
         return kTRUE;
      }
      switch ((EKeySym)keysym) {   // ignore these keys
         case kKey_Shift:
         case kKey_Control:
         case kKey_Meta:
         case kKey_Alt:
         case kKey_CapsLock:
         case kKey_NumLock:
         case kKey_ScrollLock:
            return kTRUE;
         default:
            break;
      }
      if (event->fState & kKeyControlMask) {   // Cntrl key modifier pressed
         switch ((EKeySym)keysym & ~0x20) {   // treat upper and lower the same
            case kKey_B:
               fMenuFile->Activated(kBrowse);
               return kTRUE;
            case kKey_O:
               fMenuFile->Activated(kOpenFile);
               return kTRUE;
            case kKey_E:
               fMenuFile->Activated(kNewEditor);
               return kTRUE;
            case kKey_C:
               fMenuFile->Activated(kNewCanvas);
               return kTRUE;
            case kKey_H:
               fMenuFile->Activated(kNewHtml);
               return kTRUE;
            case kKey_N:
               fMenuFile->Activated(kClone);
               return kTRUE;
            case kKey_T:
               fMenuFile->Activated(kCloseTab);
               return kTRUE;
            case kKey_W:
               fMenuFile->Activated(kCloseWindow);
               return kTRUE;
            case kKey_Q:
               fMenuFile->Activated(kQuitRoot);
               return kTRUE;
            default:
               break;
         }
      }
   }
   return TGMainFrame::HandleKey(event);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle menu entries events.

void TRootBrowser::HandleMenu(Int_t id)
{
   TRootHelpDialog *hd;
   TString cmd;
   static Int_t eNr = 1;
   TGPopupMenu *sender = (TGPopupMenu *)gTQSender;
   if (sender != fMenuFile)
      return;
   switch (id) {
      case kBrowse:
         new TBrowser();
         break;
      case kOpenFile:
         {
            Bool_t newfile = kFALSE;
            static TString dir(".");
            TGFileInfo fi;
            fi.fFileTypes = gOpenFileTypes;
            fi.SetIniDir(dir);
            new TGFileDialog(gClient->GetDefaultRoot(), this,
                             kFDOpen,&fi);
            dir = fi.fIniDir;
            if (fi.fMultipleSelection && fi.fFileNamesList) {
               TObjString *el;
               TIter next(fi.fFileNamesList);
               while ((el = (TObjString *) next())) {
                  gROOT->ProcessLine(Form("new TFile(\"%s\");",
                                     gSystem->UnixPathName(el->GetString())));
               }
               newfile = kTRUE;
            }
            else if (fi.fFilename) {
               gROOT->ProcessLine(Form("new TFile(\"%s\");",
                                  gSystem->UnixPathName(fi.fFilename)));
               newfile = kTRUE;
            }
            if (fActBrowser && newfile) {
               TGFileBrowser *fb = dynamic_cast<TGFileBrowser *>(fActBrowser);
               if (fb) fb->Selected(0);
            }
         }
         break;
                  // Handle Help menu items...
      case kHelpAbout:
         {
#ifdef R__UNIX
            TString rootx = TROOT::GetBinDir() + "/root -a &";
            gSystem->Exec(rootx);
#else
#ifdef WIN32
            new TWin32SplashThread(kTRUE);
#else
            char str[32];
            sprintf(str, "About ROOT %s...", gROOT->GetVersion());
            hd = new TRootHelpDialog(this, str, 600, 400);
            hd->SetText(gHelpAbout);
            hd->Popup();
#endif
#endif
         }
         break;
      case kHelpOnCanvas:
         hd = new TRootHelpDialog(this, "Help on Canvas...", 600, 400);
         hd->SetText(gHelpCanvas);
         hd->Popup();
         break;
      case kHelpOnMenus:
         hd = new TRootHelpDialog(this, "Help on Menus...", 600, 400);
         hd->SetText(gHelpPullDownMenus);
         hd->Popup();
         break;
      case kHelpOnGraphicsEd:
         hd = new TRootHelpDialog(this, "Help on Graphics Editor...", 600, 400);
         hd->SetText(gHelpGraphicsEditor);
         hd->Popup();
         break;
      case kHelpOnBrowser:
         hd = new TRootHelpDialog(this, "Help on Browser...", 600, 400);
         hd->SetText(gHelpBrowser);
         hd->Popup();
         break;
      case kHelpOnObjects:
         hd = new TRootHelpDialog(this, "Help on Objects...", 600, 400);
         hd->SetText(gHelpObjects);
         hd->Popup();
         break;
      case kHelpOnPS:
         hd = new TRootHelpDialog(this, "Help on PostScript...", 600, 400);
         hd->SetText(gHelpPostscript);
         hd->Popup();
         break;
      case kHelpOnRemote:
         hd = new TRootHelpDialog(this, "Help on Browser...", 600, 400);
         hd->SetText(gHelpRemote);
         hd->Popup();
         break;
      case kClone:
         CloneBrowser();
         break;
      case kNewEditor:
         cmd.Form("new TGTextEditor((const char *)0, gClient->GetRoot())");
         ++eNr;
         ExecPlugin(Form("Editor %d", eNr), "", cmd.Data(), 1);
         break;
      case kNewCanvas:
         if (IsWebGUI())
            gROOT->ProcessLine("new TCanvas()");
         else
            ExecPlugin("", "", "new TCanvas()", 1);
         break;
      case kNewHtml:
         cmd.Form("new TGHtmlBrowser(\"%s\", gClient->GetRoot())",
                  gEnv->GetValue("Browser.StartUrl", "http://root.cern.ch"));
         ExecPlugin("HTML", "", cmd.Data(), 1);
         break;
      case kExecPluginMacro:
         {
            static TString dir(".");
            TGFileInfo fi;
            fi.fFileTypes = gPluginFileTypes;
            fi.SetIniDir(dir);
            new TGFileDialog(gClient->GetDefaultRoot(), this,
                             kFDOpen,&fi);
            dir = fi.fIniDir;
            if (fi.fFilename) {
               ExecPlugin(0, fi.fFilename, 0, kRight);
            }
         }
         break;
      case kExecPluginCmd:
         {
            char command[1024];
            strlcpy(command, "new TGLSAViewer(gClient->GetRoot(), 0);",
                    sizeof(command));
            new TGInputDialog(gClient->GetRoot(), this,
                              "Enter plugin command line:",
                              command, command);
            if (strcmp(command, "")) {
               ExecPlugin("User", 0, command, kRight);
            }
         }
         break;
      case kCloseTab:
         CloseTab(fTabRight->GetCurrent());
         break;
      case kCloseWindow:
         CloseWindow();
         break;
      case kQuitRoot:
         CloseWindow();
         gApplication->Terminate(0);
         break;
      default:
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize default plugins. Could be also of the form:
/// StartEmbedding(0);
/// TPluginHandler *ph;
/// ph = gROOT->GetPluginManager()->FindHandler("TGClassBrowser");
/// if (ph && ph->LoadPlugin() != -1) {
///    ph->ExecPlugin(3, gClient->GetRoot(), 200, 500);
/// }
/// StopEmbedding();

void TRootBrowser::InitPlugins(Option_t *opt)
{
   TString cmd;

   if ((opt == 0) || (!opt[0]))
      return;
   // --- Left vertical area

   // File Browser plugin
   if (strchr(opt, 'F')) {
      cmd.Form("new TGFileBrowser(gClient->GetRoot(), (TBrowser *)0x%lx, 200, 500);", (ULong_t)fBrowser);
      ExecPlugin("Files", 0, cmd.Data(), 0);
      ++fNbInitPlugins;
   }

   // --- Right main area

   Int_t i, len = strlen(opt);
   for (i=0; i<len; ++i) {
      // Editor plugin...
      if (opt[i] == 'E') {
         cmd.Form("new TGTextEditor((const char *)0, gClient->GetRoot());");
         ExecPlugin("Editor 1", 0, cmd.Data(), 1);
         ++fNbInitPlugins;
      }

      // HTML plugin...
      if (opt[i] == 'H') {
         if (gSystem->Load("libGuiHtml") >= 0) {
            cmd.Form("new TGHtmlBrowser(\"%s\", gClient->GetRoot());",
                     gEnv->GetValue("Browser.StartUrl",
                     "http://root.cern.ch/root/html/ClassIndex.html"));
            ExecPlugin("HTML", 0, cmd.Data(), 1);
            ++fNbInitPlugins;
         }
      }

      // Canvas plugin...
      if ((opt[i] == 'C') && !IsWebGUI()) {
         cmd.Form("new TCanvas();");
         ExecPlugin("c1", 0, cmd.Data(), 1);
         ++fNbInitPlugins;
      }

      // GLViewer plugin...
      if (opt[i] == 'G') {
         cmd.Form("new TGLSAViewer(gClient->GetRoot(), 0);");
         ExecPlugin("OpenGL", 0, cmd.Data(), 1);
         ++fNbInitPlugins;
      }

      // PROOF plugin...
      if (opt[i] == 'P') {
         cmd.Form("new TSessionViewer();");
         ExecPlugin("PROOF", 0, cmd.Data(), 1);
         ++fNbInitPlugins;
      }
   }
   // --- Right bottom area

   // Command plugin...
   if (strchr(opt, 'I')) {
      cmd.Form("new TGCommandPlugin(gClient->GetRoot(), 700, 300);");
      ExecPlugin("Command", 0, cmd.Data(), 2);
      ++fNbInitPlugins;
   }

   // --- Select first tab everywhere
   SetTab(0, 0);
   SetTab(1, 0);
   SetTab(2, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the GUI factory is set to use the Web GUI.

Bool_t TRootBrowser::IsWebGUI()
{
   TString factory = gEnv->GetValue("Gui.Factory", "native");
   return (factory.Contains("web", TString::kIgnoreCase));
}

////////////////////////////////////////////////////////////////////////////////
/// Really delete the browser and the this GUI.

void TRootBrowser::ReallyDelete()
{
   if (fBrowser) fBrowser->SetBrowserImp(0);
   delete this;
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively remove object from browser.

void TRootBrowser::RecursiveRemove(TObject *obj)
{
   if (fActBrowser)
      fActBrowser->RecursiveRemove(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively reparent TGPopupMenu to gClient->GetDefaultRoot().

void TRootBrowser::RecursiveReparent(TGPopupMenu *popup)
{
   TGMenuEntry *entry = 0;
   TIter next(popup->GetListOfEntries());
   while ((entry = (TGMenuEntry *)next())) {
      if (entry->GetPopup()) {
         RecursiveReparent(entry->GetPopup());
      }
   }
   popup->ReparentWindow(gClient->GetDefaultRoot());
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh the actual browser contents.

void TRootBrowser::Refresh(Bool_t force)
{
   if (fActBrowser)
      fActBrowser->Refresh(force);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove tab element "subpos" from tab "pos".

void TRootBrowser::RemoveTab(Int_t pos, Int_t subpos)
{
   TGTab *edit = 0;
   switch (pos) {
      case kLeft: // left
         edit = fTabLeft;
         break;
      case kRight: // right
         edit = fTabRight;
         fMenuFrame->HideFrame(fActMenuBar);
         fMenuFrame->GetList()->Remove(fActMenuBar);
         fActMenuBar = 0;
         break;
      case kBottom: // bottom
         edit = fTabBottom;
         break;
   }
   if (!edit || !edit->GetTabTab(subpos))
      return;
   const char *tabName = edit->GetTabTab(subpos)->GetString();
   TObject *obj = 0;
   if ((obj = fPlugins.FindObject(tabName))) {
      fPlugins.Remove(obj);
   }
   TGFrameElement *el = 0;
   if (edit->GetTabContainer(subpos))
      el = (TGFrameElement *)edit->GetTabContainer(subpos)->GetList()->First();
   if (el && el->fFrame) {
      el->fFrame->Disconnect("ProcessedConfigure(Event_t*)");
      el->fFrame->SetFrameElement(0);
      if (el->fFrame->InheritsFrom("TGMainFrame")) {
         Bool_t sleep = (el->fFrame->InheritsFrom("TRootCanvas")) ? kTRUE : kFALSE;
         ((TGMainFrame *)el->fFrame)->CloseWindow();
         if (sleep)
            gSystem->Sleep(150);
         gSystem->ProcessEvents();
      }
      else
         delete el->fFrame;
      el->fFrame = 0;
      if (el->fLayout && (el->fLayout != fgDefaultHints) &&
         (el->fLayout->References() > 0)) {
         el->fLayout->RemoveReference();
         if (!el->fLayout->References()) {
            delete el->fLayout;
         }
      }
      edit->GetTabContainer(subpos)->GetList()->Remove(el);
      delete el;
   }
   fNbTab[pos]--;
   edit->RemoveTab(subpos);
   SwitchMenus(edit->GetTabContainer(edit->GetCurrent()));
}

////////////////////////////////////////////////////////////////////////////////
/// Switch to Tab "subpos" in TGTab "pos".

void TRootBrowser::SetTab(Int_t pos, Int_t subpos)
{
   TGTab *tab = GetTab(pos);
   if (subpos == -1)
      subpos = fCrTab[pos];

   if (tab && tab->SetTab(subpos, kFALSE)) { // Block signal emit
      if (pos == kRight)
         SwitchMenus(tab->GetTabContainer(subpos));
      tab->Layout();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set text "title" of Tab "subpos" in TGTab "pos".

void TRootBrowser::SetTabTitle(const char *title, Int_t pos, Int_t subpos)
{
   TBrowserPlugin *p = 0;
   TGTab *edit = GetTab(pos);
   if (!edit) return;
   if (subpos == -1)
      subpos = fCrTab[pos];

   TGTabElement *el = edit->GetTabTab(subpos);
   if (el) {
      el->SetText(new TGString(title));
      edit->Layout();
      if ((p = (TBrowserPlugin *)fPlugins.FindObject(title)))
         p->SetName(title);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set text in culumn col in status bar.

void TRootBrowser::SetStatusText(const char* txt, Int_t col)
{
   fStatusBar->SetText(txt, col);
}

////////////////////////////////////////////////////////////////////////////////
/// Show the selected frame's menu and hide previous one.

void TRootBrowser::ShowMenu(TGCompositeFrame *menu)
{
   TGFrameElement *el = 0;
   // temporary solution until I find a proper way to handle
   // these bloody menus...
   fBindList->Delete();
   TIter nextm(fMenuBar->GetList());
   while ((el = (TGFrameElement *) nextm())) {
      TGMenuTitle *t = (TGMenuTitle *) el->fFrame;
      Int_t code = t->GetHotKeyCode();
      BindKey(fMenuBar, code, kKeyMod1Mask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyShiftMask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyLockMask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyMod2Mask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
   }
   fMenuFrame->HideFrame(fActMenuBar);
   fMenuFrame->ShowFrame(menu);
   menu->Layout();
   fMenuFrame->Layout();
   fActMenuBar = menu;
}

////////////////////////////////////////////////////////////////////////////////
/// Start embedding external frame in the tab "pos" and tab element "subpos".

void TRootBrowser::StartEmbedding(Int_t pos, Int_t subpos)
{
   fEditTab = GetTab(pos);
   if (!fEditTab) return;
   fEditPos = pos;
   fEditSubPos = subpos;

   if (fEditFrame == 0) {
      if (subpos == -1) {
         fCrTab[pos] = fNbTab[pos]++;
         fEditFrame  = fEditTab->AddTab(Form("Tab %d",fNbTab[pos]));
         fEditSubPos = fEditTab->GetNumberOfTabs()-1;
         fEditFrame->MapWindow();
         TGTabElement *tabel = fEditTab->GetTabTab(fEditSubPos);
         if(tabel) {
            tabel->MapWindow();
            if (fShowCloseTab && (pos == 1))
               tabel->ShowClose();
         }
         fEditTab->SetTab(fEditTab->GetNumberOfTabs()-1);
         fEditTab->Layout();
      }
      else {
         fCrTab[pos] = subpos;
         fEditFrame = fEditTab->GetTabContainer(subpos);
         fEditTab->SetTab(subpos);
      }
      if (fEditFrame) fEditFrame->SetEditable();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stop embedding external frame in the current editable frame.

void TRootBrowser::StopEmbedding(const char *name, TGLayoutHints *layout)
{
   if (fEditFrame != 0) {
      fEditFrame->SetEditable(kFALSE);
      TGFrameElement *el = (TGFrameElement*) fEditFrame->GetList()->First();
      if (el && el->fFrame) {
         // let be notified when the inside frame gets resized, and tell its
         // container to recompute its layout
         el->fFrame->Connect("ProcessedConfigure(Event_t*)", "TGCompositeFrame",
                             fEditFrame, "Layout()");
      }
      if (layout) {
         el = (TGFrameElement*) fEditFrame->GetList()->Last();
         // !!!! MT what to do with the old layout? Leak it for now ...
         if (el) el->fLayout = layout;
      }
      fEditFrame->Layout();
      if (fEditTab == fTabRight)
         SwitchMenus(fEditFrame);
   }
   if (name && strlen(name)) {
      SetTabTitle(name, fEditPos, fEditSubPos);
   }
   if (fEditTab) fEditTab->Selected(fEditSubPos);
   fEditFrame = fEditTab = 0;
   fEditPos = fEditSubPos = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Move the menu from original frame to our TGMenuFrame, or display the
/// menu associated to the current tab.

void TRootBrowser::SwitchMenus(TGCompositeFrame  *from)
{
   if (from == 0)
      return;
   TGFrameElement *fe = (TGFrameElement *)from->GetList()->First();
   if (!fe) {
      if (fActMenuBar != fMenuBar)
         ShowMenu(fMenuBar);
      return;
   }
   TGCompositeFrame *embed = (TGCompositeFrame *)fe->fFrame;
   TGFrameElement *el = 0;
   if (embed && embed->GetList()) {
      TIter next(embed->GetList());
      while ((el = (TGFrameElement *)next())) {
         if (el->fFrame->InheritsFrom("TGMenuBar")) {
            TGMenuBar *menu = (TGMenuBar *)el->fFrame;
            if (fActMenuBar == menu)
               return;
            TGFrameElement *nw;
            TIter nel(fMenuFrame->GetList());
            while ((nw = (TGFrameElement *) nel())) {
               if (nw->fFrame == menu) {
                  ShowMenu(menu);
                  return;
               }
            }
            ((TGCompositeFrame *)menu->GetParent())->HideFrame(menu);
            ((TGCompositeFrame *)menu->GetParent())->SetCleanup(kNoCleanup);
            menu->ReparentWindow(fMenuFrame);
            fMenuFrame->AddFrame(menu, fLH2);
            TGFrameElement *mel;
            TIter mnext(menu->GetList());
            while ((mel = (TGFrameElement *) mnext())) {
               TGMenuTitle *t = (TGMenuTitle *) mel->fFrame;
               TGPopupMenu *popup = menu->GetPopup(t->GetName());
               if (popup) {
                  RecursiveReparent(popup);
                  if (popup->GetEntry("Close Canvas")) {
                     TGMenuEntry *exit = popup->GetEntry("Close Canvas");
                     popup->HideEntry(exit->GetEntryId());
                  }
                  if (popup->GetEntry("Close Viewer")) {
                     TGMenuEntry *exit = popup->GetEntry("Close Viewer");
                     popup->HideEntry(exit->GetEntryId());
                  }
                  if (popup->GetEntry("Quit ROOT")) {
                     TGMenuEntry *exit = popup->GetEntry("Quit ROOT");
                     popup->HideEntry(exit->GetEntryId());
                  }
                  if (popup->GetEntry("Exit")) {
                     TGMenuEntry *exit = popup->GetEntry("Exit");
                     popup->HideEntry(exit->GetEntryId());
                  }
               }
            }
            ShowMenu(menu);
            return;
         }
      }
   }
   if (fActMenuBar != fMenuBar)
      ShowMenu(fMenuBar);
}

////////////////////////////////////////////////////////////////////////////////
/// Emits signal when double clicking on icon.

void TRootBrowser::DoubleClicked(TObject *obj)
{
   Emit("DoubleClicked(TObject*)", (Long_t)obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Emits signal when double clicking on icon.

void TRootBrowser::Checked(TObject *obj, Bool_t checked)
{
   Long_t args[2];

   args[0] = (Long_t)obj;
   args[1] = checked;

   Emit("Checked(TObject*,Bool_t)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emits signal "ExecuteDefaultAction(TObject*)".

void TRootBrowser::ExecuteDefaultAction(TObject *obj)
{
   Emit("ExecuteDefaultAction(TObject*)", (Long_t)obj);
}


////////////////////////////////////////////////////////////////////////////////
/// static contructor returning TBrowserImp,
/// as needed by the plugin mechanism.

TBrowserImp *TRootBrowser::NewBrowser(TBrowser *b, const char *title,
                                      UInt_t width, UInt_t height,
                                      Option_t *opt)
{
   TRootBrowser *browser = new TRootBrowser(b, title, width, height, opt);
   return (TBrowserImp *)browser;
}

////////////////////////////////////////////////////////////////////////////////
/// static contructor returning TBrowserImp,
/// as needed by the plugin mechanism.

TBrowserImp *TRootBrowser::NewBrowser(TBrowser *b, const char *title, Int_t x,
                                      Int_t y, UInt_t width, UInt_t height,
                                      Option_t *opt)
{
   TRootBrowser *browser = new TRootBrowser(b, title, x, y, width, height, opt);
   return (TBrowserImp *)browser;
}
