// @(#)root/gui:$Id$
// Author: Bertrand Bellenot   26/09/2007

#include "TROOT.h"
#include "TSystem.h"
#include "TRint.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TGLabel.h"
#include "TGFrame.h"
#include "TGLayout.h"
#include "TGComboBox.h"
#include "TGTextView.h"
#include "TGTextEntry.h"
#include "TInterpreter.h"
#include "Getline.h"

#include "TGCommandPlugin.h"

//_____________________________________________________________________________
//
// TGCommandPlugin
//
// Class used to redirect command line input/output.
//_____________________________________________________________________________

ClassImp(TGCommandPlugin);

////////////////////////////////////////////////////////////////////////////////
/// TGCommandPlugin Constructor.

TGCommandPlugin::TGCommandPlugin(const TGWindow *p, UInt_t w, UInt_t h) :
      TGMainFrame(p, w, h)
{
   SetCleanup(kDeepCleanup);
   fHf = new TGHorizontalFrame(this, 100, 20);
   fComboCmd   = new TGComboBox(fHf, "", 1);
   fCommand    = fComboCmd->GetTextEntry();
   fCommandBuf = fCommand->GetBuffer();
   fComboCmd->Resize(200, fCommand->GetDefaultHeight());
   fHf->AddFrame(fComboCmd, new TGLayoutHints(kLHintsCenterY |
                 kLHintsRight | kLHintsExpandX, 5, 5, 1, 1));
   fHf->AddFrame(fLabel = new TGLabel(fHf, "Command (local):"),
                 new TGLayoutHints(kLHintsCenterY | kLHintsRight,
                 5, 5, 1, 1));
   AddFrame(fHf, new TGLayoutHints(kLHintsLeft | kLHintsTop |
            kLHintsExpandX, 3, 3, 3, 3));
   fCommand->Connect("ReturnPressed()", "TGCommandPlugin", this,
                     "HandleCommand()");
   fStatus = new TGTextView(this, 10, 100, 1);
   if (gClient->GetStyle() < 2) {
      Pixel_t pxl;
      gClient->GetColorByName("#a0a0a0", pxl);
      fStatus->SetSelectBack(pxl);
      fStatus->SetSelectFore(TGFrame::GetWhitePixel());
   }
   AddFrame(fStatus, new TGLayoutHints(kLHintsLeft | kLHintsTop |
            kLHintsExpandX | kLHintsExpandY, 3, 3, 3, 3));
   fPid = gSystem->GetPid();
   TString defhist(Form("%s/.root_hist", gSystem->UnixPathName(
                        gSystem->HomeDirectory())));
   FILE *lunin = fopen(defhist.Data(), "rt");
   if (lunin) {
      ULong_t linecount = 0;
      char histline[256];
      rewind(lunin);
      while (fgets(histline, 256, lunin))
         ++linecount;
      rewind(lunin);
      if (linecount > 500) {
         linecount -= 500;
         while(--linecount > 0)
            if (!fgets(histline, 256, lunin))
               break;
      }
      linecount = 0;
      while (fgets(histline, 256, lunin)) {
         histline[strlen(histline)-1] = 0; // remove trailing "\n"
         fComboCmd->InsertEntry(histline, 0, -1);
         // limit the history size to 500 lines
         if (++linecount > 500)
            break;
      }
      fclose(lunin);
   }
   fTimer = new TTimer(this, 1000);
   fTimer->Reset();
   fTimer->TurnOn();
   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGCommandPlugin::~TGCommandPlugin()
{
   TString pathtmp = TString::Format("%s/command.%d.log",
                                     gSystem->TempDirectory(), fPid);
   gSystem->Unlink(pathtmp);
   fCommand->Disconnect("ReturnPressed()");
   delete fTimer;
   fTimer = 0;
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if actual ROOT session is a remote one or a local one.

void TGCommandPlugin::CheckRemote(const char * /*str*/)
{
   Pixel_t pxl;
   TApplication *app = gROOT->GetApplication();
   if (!app->InheritsFrom("TRint"))
      return;
   TString sPrompt = ((TRint*)app)->GetPrompt();
   Int_t end = sPrompt.Index(":root [", 0);
   if (end > 0 && end != kNPOS) {
      // remote session
      sPrompt.Remove(end);
      gClient->GetColorByName("#ff0000", pxl);
      fLabel->SetTextColor(pxl);
      fLabel->SetText(Form("Command (%s):", sPrompt.Data()));
   }
   else {
      // local session
      gClient->GetColorByName("#000000", pxl);
      fLabel->SetTextColor(pxl);
      fLabel->SetText("Command (local):");
   }
   fHf->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle command line from the "command" combo box.

void TGCommandPlugin::HandleCommand()
{
   const char *string = fCommandBuf->GetString();
   if (strlen(string) > 0) {
      // form temporary file path
      TString sPrompt = "root []";
      TString pathtmp = TString::Format("%s/command.%d.log",
                                        gSystem->TempDirectory(), fPid);
      TApplication *app = gROOT->GetApplication();
      if (app->InheritsFrom("TRint"))
         sPrompt = ((TRint*)gROOT->GetApplication())->GetPrompt();
      FILE *lunout = fopen(pathtmp.Data(), "a+t");
      if (lunout) {
         fputs(Form("%s%s\n",sPrompt.Data(), string), lunout);
         fclose(lunout);
      }
      gSystem->RedirectOutput(pathtmp.Data(), "a");
      gApplication->SetBit(TApplication::kProcessRemotely);
      gROOT->ProcessLine(string);
      fComboCmd->InsertEntry(string, 0, -1);
      if (app->InheritsFrom("TRint"))
         Gl_histadd((char *)string);
      gSystem->RedirectOutput(0);
      fStatus->LoadFile(pathtmp.Data());
      fStatus->ShowBottom();
      CheckRemote(string);
      fCommand->Clear();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle timer event.

Bool_t TGCommandPlugin::HandleTimer(TTimer *t)
{
   if ((fTimer == 0) || (t != fTimer)) return kTRUE;
   CheckRemote("");
   return kTRUE;
}
