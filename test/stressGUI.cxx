// @(#)root/test:$name:  $:$id: stressGUI.cxx,v 1.0 exp $
// Author: Bertrand Bellenot

//
//    ROOT GUI test suite.
//
// The suite of programs below tests many elements of the ROOT GUI classes
//
// The test can only be run as a standalone program.
// To build and run it:
//
//    make stressGUI
//    stressGUI
//
// To get a short help:
//    stressGUI -help
//

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <snprintf.h>
#include <TString.h>
#include <TROOT.h>
#include <TClass.h>
#include <TEnv.h>
#include <TError.h>
#include <TBenchmark.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TDatime.h>
#include <TFile.h>
#include <TObjArray.h>

#include <TGFrame.h>
#include <TImage.h>
#include <TCanvas.h>

#include <TMD5.h>
#include <TG3DLine.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGColorDialog.h>
#include <TGComboBox.h>
#include <TGLabel.h>
#include <TGListBox.h>
#include <TGListTree.h>
#include <TGMenu.h>
#include <TGMsgBox.h>
#include <TGNumberEntry.h>
#include <TGProgressBar.h>
#include <TGResourcePool.h>
#include <TGShutter.h>
#include <TGSimpleTable.h>
#include <TGTextEdit.h>
#include <TRootCanvas.h>
#include <TGTab.h>
#include <TGPack.h>
#include <TGColorDialog.h>
#include <TGFontDialog.h>
#include <TGTextEditDialogs.h>
#include <TGTableLayout.h>
#include <TGMdi.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>
#include <TGTripleSlider.h>
#include <TBrowser.h>
#include <TGPasswdDialog.h>
#include <TGImageMap.h>
#include <TASPaletteEditor.h>
#include <TControlBar.h>
#include <TControlBarImp.h>
#include <TGSpeedo.h>
#include <TGShapedFrame.h>
#include <TGSplitFrame.h>
#include <TGTextEditor.h>
#include <TRootHelpDialog.h>
#include <TGHtmlBrowser.h>
#include <HelpText.h>
#include <TSystemDirectory.h>
#include <TInterpreter.h>
#include <TStopwatch.h>
#include <TVirtualX.h>

#include <TRecorder.h>

void     stressGUI();
void     ProcessFrame(TGFrame *f, const char *title);

// Tests functions.
void     testLayout();
void     testTextAlign();
void     testGroupState();
void     testLabels();
void     testSplitButton();
void     testTextEntries();
void     testListTree();
void     testShutter();
void     testProgressBar();
void     testNumberEntry();
void     testEditor();
void     testCanvas();
void     testColorDlg();
void     testFontDlg();
void     testSearchDlg();
void     testTableLayout();
void     testPack();
void     testSliders();
void     testBrowsers();
void     testSplitFrame();
void     testControlBars();
void     testHelpDialog();
void     testPaletteEditor();
void     testHtmlBrowser();

void     run_tutorials();
void     guitest_playback();
void     dnd_playback();
void     mditest_playback();
void     fitpanel_playback();
void     graph_edit_playback();

// Global variables.
RedirectHandle_t gRH;
Int_t    gTestNum = 0;
Bool_t   gOptionRef  = kFALSE;
Bool_t   gOptionKeep = kFALSE;
Bool_t   gOptionFull = kFALSE;
char     outfile[80];
char     gLine[80];
Int_t    sizes[100];
TString  gTmpfilename;
TString  gRootSys;

FILE    *sgref = 0;

////////////////////////////////////////////////////////////////////////////////
/// Application main entry point.

int main(int argc, char *argv[])
{
   // use $ROOTSYS/etc/system.rootrc default values
   gEnv->ReadFile(TString::Format("%s/etc/system.rootrc",
                  gSystem->Getenv("ROOTSYS")), kEnvAll);
   gOptionRef  = kFALSE;
   gOptionKeep = kFALSE;
   gOptionFull = kFALSE;
   gTmpfilename = "stress-gui";
   FILE *f = gSystem->TempFileName(gTmpfilename);
   fclose(f);
   for (int i = 0; i < argc; i++) {
      if (!strcmp(argv[i], "-ref")) gOptionRef = kTRUE;
      if (!strcmp(argv[i], "-keep")) gOptionKeep = kTRUE;
      if (!strcmp(argv[i], "-full")) gOptionFull = kTRUE;
      if (!strcmp(argv[i], "-help") || !strcmp(argv[i], "-?")) {
         printf("Usage: stressGUI [-ref] [-keep] [-full] [-help] [-?] \n");
         printf("Options:\n");
         printf("\n");
         printf("  -ref: Generate the reference output file \"stressGUI.ref\"\n");
         printf("\n");
         printf("  -keep: Keep the png files even for passed tests\n");
         printf("        (by default the png files are deleted)\n");
         printf("\n");
         printf("  -full: Full test: replay also recorder sessions\n");
         printf("        (guitest, drag and drop, fitpanel, ...)\n");
         printf("\n");
         printf("  -help, -?: Print usage and exit\n");
         return 0;
      }
   }
   TApplication theApp("App", &argc, argv);
   gBenchmark = new TBenchmark();
   stressGUI();
   theApp.Terminate();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Run all stress GUI tests.

void stressGUI()
{
   if (gOptionRef) {
      sgref = fopen("stressGUI.ref", "wt");
   }
   else {
      // Read the reference file "stressGUI.ref"
      sgref = fopen("stressGUI.ref", "rt");
      if (sgref == 0) {
         printf("\nReference file \"stressGUI.ref\" not found!\n");
         printf("Please generate the reference file by executing\n");
         printf("stressGUI with the -ref flag, as shown below:\n");
         printf("   stressGUI -ref\n");
         gSystem->Unlink(gTmpfilename.Data());
         exit(0);
      }
      char line[160];
      Int_t i = -1;
      while (fgets(line, 160, sgref)) {
         if ((i >= 0) && (strlen(line) > 15)) {
            sscanf(&line[8],  "%d", &sizes[i]);
         }
         i++;
      }
      fclose(sgref);
   }
   gRootSys = gSystem->UnixPathName(gSystem->Getenv("ROOTSYS"));
#ifdef WIN32
   // remove the drive letter (e.g. "C:/") from $ROOTSYS, if any
   if (gRootSys[1] == ':' && gRootSys[2] == '/')
      gRootSys.Remove(0, 2);
#endif

   gVirtualX->Warp(gClient->GetDisplayWidth()-50, gClient->GetDisplayHeight()-50,
                   gClient->GetDefaultRoot()->GetId());
   // uncomment the next few lines to avoid (forbid) any mouse interaction
//   gVirtualX->GrabPointer(gClient->GetDefaultRoot()->GetId(), kButtonPressMask |
//                          kButtonReleaseMask | kPointerMotionMask, kNone,
//                          gVirtualX->CreateCursor(kWatch), kTRUE, kFALSE);

   if (gOptionRef) {
      fprintf(sgref, "Test#     Size#\n");
   } else {
      std::cout << "**********************************************************************" <<std::endl;
      std::cout << "*  Starting  GUI - S T R E S S suite                                 *" <<std::endl;
      std::cout << "**********************************************************************" <<std::endl;
   }
   gTestNum = 0;

   gBenchmark->Start("stressGUI");

   if (!gOptionRef) {
      std::cout << "*  Running macros in $ROOTSYS/tutorials/gui - S T R E S S            *" <<std::endl;
      std::cout << "**********************************************************************" <<std::endl;
   }
   run_tutorials();
   if (!gOptionRef) {
      std::cout << "**********************************************************************" <<std::endl;
      std::cout << "*  Starting Basic GUI Widgets - S T R E S S                          *" <<std::endl;
      std::cout << "**********************************************************************" <<std::endl;
   }
   testLayout();
   testTextAlign();
   testGroupState();
   testLabels();
   testSplitButton();
   testTextEntries();
   testListTree();
   testShutter();
   testProgressBar();
   testNumberEntry();
   testTableLayout();
   if (!gOptionRef) {
      std::cout << "**********************************************************************" <<std::endl;
      std::cout << "*  Starting High Level GUI Widgets - S T R E S S                     *" <<std::endl;
      std::cout << "**********************************************************************" <<std::endl;
   }
   testPack();
   testSearchDlg();
   testFontDlg();
   testColorDlg();
   testEditor();
   testCanvas();
   testSliders();
   testBrowsers();
   testSplitFrame();
   testControlBars();
   testHelpDialog();
   testPaletteEditor();
   testHtmlBrowser();

   if (!gOptionRef) {

      if (gOptionFull) {
         std::cout << "**********************************************************************" <<std::endl;
         std::cout << "*  Starting Drag and Drop playback - S T R E S S                     *" <<std::endl;
         std::cout << "**********************************************************************" <<std::endl;
         dnd_playback();

         std::cout << "**********************************************************************" <<std::endl;
         std::cout << "*  Starting MDI test playback - S T R E S S                          *" <<std::endl;
         std::cout << "**********************************************************************" <<std::endl;
         mditest_playback();

         std::cout << "**********************************************************************" <<std::endl;
         std::cout << "*  Starting guitest recorder playback - S T R E S S                  *" <<std::endl;
         std::cout << "**********************************************************************" <<std::endl;
         guitest_playback();

         std::cout << "**********************************************************************" <<std::endl;
         std::cout << "*  Starting fit panel recorder playback - S T R E S S                *" <<std::endl;
         std::cout << "**********************************************************************" <<std::endl;
         fitpanel_playback();

         std::cout << "**********************************************************************" <<std::endl;
         std::cout << "*  Starting graphic editors recorder playback - S T R E S S          *" <<std::endl;
         std::cout << "**********************************************************************" <<std::endl;
         graph_edit_playback();
      }
      std::cout << "**********************************************************************" <<std::endl;

      gBenchmark->Stop("stressGUI");

      //Print table with results
      Bool_t UNIX = strcmp(gSystem->GetName(), "Unix") == 0;
      if (UNIX) {
         TString sp = gSystem->GetFromPipe("uname -a");
         sp.Resize(60);
         printf("*  SYS: %s\n",sp.Data());
         if (strstr(gSystem->GetBuildNode(),"Linux")) {
            sp = gSystem->GetFromPipe("lsb_release -d -s");
            printf("*  SYS: %s\n",sp.Data());
         }
         if (strstr(gSystem->GetBuildNode(),"Darwin")) {
            sp  = gSystem->GetFromPipe("sw_vers -productVersion");
            sp += " Mac OS X ";
            printf("*  SYS: %s\n",sp.Data());
         }
      } else {
         const char *os = gSystem->Getenv("OS");
         if (!os) printf("*  SYS: Windows 95\n");
         else     printf("*  SYS: %s %s \n",os,gSystem->Getenv("PROCESSOR_IDENTIFIER"));
      }

      printf("**********************************************************************\n");
      printf("*  ");
      gBenchmark->Print("stressGUI");

      Double_t ct = gBenchmark->GetCpuTime("stressGUI");  // ref: 13 s
      Double_t rt = gBenchmark->GetRealTime("stressGUI"); // ref: 300 s
      // normalize at 1000 rootmarks
      Double_t full_marks = 0.5 *((13.0/ct) + (300.0/rt));
      if (!gOptionFull)
         full_marks = 0.5 *((4.5/ct) + (35.0/rt));
      const Double_t rootmarks = 1000.0 * full_marks;

      printf("**********************************************************************\n");
      printf("*  ROOTMARKS = %6.1f   *  Root%-8s  %d/%04d\n", rootmarks, gROOT->GetVersion(),
             gROOT->GetVersionDate(), gROOT->GetVersionTime());
      printf("**********************************************************************\n");
   }
   gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer
   if (gOptionRef) {
      fclose(sgref);
   }
   gSystem->Unlink(gTmpfilename.Data());
#ifdef WIN32
   gSystem->Exec("erase /q /s TxtEdit* >nul 2>&1");
   gSystem->Exec("erase /q /s TxtView* >nul 2>&1");
#else
   gSystem->Exec("rm -f TxtEdit*");
   gSystem->Exec("rm -f TxtView*");
#endif
}

////////////////////////////////////////////////////////////////////////////////
//                                Utilities
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Return the size of the file "filename".

Int_t FileSize(const char *filename)
{
   FileStat_t fs;
   if (gSystem->GetPathInfo(filename, fs) == 0)
      return (Int_t)fs.fSize;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Verify the file size.

Bool_t VerifySize(const char *filename, const char *title)
{
   Int_t  ftol = 0;
   Bool_t success = kFALSE;
   Int_t fsize = FileSize(filename);
   // results depends on the files in $ROOTSYS/test
   if (strstr(title, "Browser")) ftol = 350;
   // not relevant as the CPU load and the background may vary...
   if (strstr(title, "CPUMeter.C")) ftol = 50000;
   if (strstr(title, "games.C")) ftol = 150;
   if (strstr(title, "ntupleTableTest.C")) ftol = 100;

   if (!gOptionRef) {
      if ((fsize < sizes[gTestNum] - ftol) || (fsize > sizes[gTestNum] + ftol))
         success = kFALSE;
      else
         success = kTRUE;

      snprintf(gLine,80,"Test %2d: %s", gTestNum, title);
      const Int_t nch = strlen(gLine);
      if (success) {
         std::cout << gLine;
         for (Int_t i = nch; i < 67; i++) std::cout << ".";
         std::cout << " OK" << std::endl;
      } else {
         std::cout << gLine;
         for (Int_t i = nch; i < 63; i++) std::cout << ".";
         std::cout << " FAILED" << std::endl;
         std::cout << "         File Size = "  << fsize << std::endl;
         std::cout << "          Ref Size = "  << sizes[gTestNum] << std::endl;
      }
   } else {
      fprintf(sgref, "%5d%10d\n", gTestNum, fsize);
      success = kTRUE;
   }
   if (!gOptionKeep && success) gSystem->Unlink(filename);
   return success;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a capture of frame f in a png file.

void ProcessFrame(TGFrame *f, const char *title)
{
   gSystem->ProcessEvents();
   gSystem->Sleep(250);
   gSystem->ProcessEvents();
   gErrorIgnoreLevel = 9999;

   if (gOptionRef)
      snprintf(outfile,80, "sgui_%02d_ref.png", gTestNum);
   else
      snprintf(outfile,80, "sgui_%02d.png", gTestNum);

   TImage *img = TImage::Create();
   f->RaiseWindow();
   img->FromWindow(f->GetId());
   img->WriteImage(outfile);

   if (!gOptionRef) {
      if (!strstr(title, "Pack Frames") &&
          !strstr(title, "HTML Browser")) {
         gSystem->RedirectOutput(gTmpfilename.Data(), "w", &gRH);
         ((TGMainFrame *)f)->SaveSource(Form("sgui_%02d.C", gTestNum));
         gSystem->Unlink(Form("sgui_%02d.C", gTestNum));
         gSystem->RedirectOutput(0, 0, &gRH);
      }
   }
   VerifySize(outfile, title);

   delete img;
   gErrorIgnoreLevel = 0;
   gTestNum++;
}

////////////////////////////////////////////////////////////////////////////////
/// Verify the size of a png file generated from a root macro.

void ProcessMacro(const char *macro, const char *title)
{
   Int_t   nbpass = 1, npass = 0;
   TString capture = macro;
   capture.ReplaceAll(".C", "_0.png");
   if (strstr(macro, "games.C")) nbpass = 3;
   if (strstr(macro, "galaxy_image.C")) nbpass = 2;

   while (npass < nbpass) {
      ++npass;
      capture.ReplaceAll(TString::Format("_%d.png", npass-1),
                         TString::Format("_%d.png", npass));
      VerifySize(capture.Data(), title);
      gTestNum++;
   }
}

////////////////////////////////////////////////////////////////////////////////

void CloseMainframes()
{
   TClass* clGMainFrame = TClass::GetClass("TGMainFrame");
   TGWindow* win = 0;
   TIter iWin(gClient->GetListOfWindows());
   while ((win = (TGWindow*)iWin())) {
      const TObject* winGetParent = win->GetParent();
      Bool_t winIsMapped = kFALSE;
      if (winGetParent == gClient->GetDefaultRoot())
         winIsMapped = win->IsMapped();
      if (winIsMapped && win->InheritsFrom(clGMainFrame)) {
         ((TGMainFrame *)win)->CloseWindow();
         gSystem->Sleep(100);
      }
      gSystem->ProcessEvents();
   }
   gSystem->Sleep(100);
}

////////////////////////////////////////////////////////////////////////////////
//                            GUI Test code
////////////////////////////////////////////////////////////////////////////////

class ButtonLayoutWindow : public TGMainFrame {

private:
   TGTextButton *test, *draw, *help, *ok, *cancel, *exit;

public:
   ButtonLayoutWindow(const TGWindow *p, UInt_t w, UInt_t h);

};

////////////////////////////////////////////////////////////////////////////////
/// Create a container frames containing buttons

ButtonLayoutWindow::ButtonLayoutWindow(const TGWindow *p, UInt_t w, UInt_t h) :
   TGMainFrame(p, w, h)
{
   SetCleanup(kDeepCleanup);
   // one button is resized up to the parent width. Note! this width should be fixed!
   TGVerticalFrame *hframe1 = new TGVerticalFrame(this, 170, 50, kFixedWidth);
   test = new TGTextButton(hframe1, "&Test ");
   // to take whole space we need to use kLHintsExpandX layout hints
   hframe1->AddFrame(test, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2,0,2,2));
   AddFrame(hframe1, new TGLayoutHints(kLHintsCenterX, 2, 2, 5, 1));

   // two buttons are resized up to the parent width. Note! this width should be fixed!
   TGCompositeFrame *cframe1 = new TGCompositeFrame(this, 170, 20, kHorizontalFrame | kFixedWidth);
   draw = new TGTextButton(cframe1, "&Draw");
   // to share whole parent space we need to use kLHintsExpandX layout hints
   cframe1->AddFrame(draw, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2,2,2,2));

   // button background will be set to yellow
   ULong_t yellow;
   gClient->GetColorByName("yellow", yellow);
   help = new TGTextButton(cframe1, "&Help");
   help->ChangeBackground(yellow);
   cframe1->AddFrame(help, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2,2,2,2));
   AddFrame(cframe1, new TGLayoutHints(kLHintsCenterX, 2, 2, 5, 1));

   // three buttons are resized up to the parent width. Note! this width should be fixed!
   TGCompositeFrame *cframe2 = new TGCompositeFrame(this, 170, 20, kHorizontalFrame | kFixedWidth);
   ok = new TGTextButton(cframe2, "OK");
   ok->SetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
   ok->SetEnabled(kFALSE);
   // to share whole parent space we need to use kLHintsExpandX layout hints
   cframe2->AddFrame(ok, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,3,2,2,2));

   TGGC myGC = *gClient->GetResourcePool()->GetFrameGC();
   //TGFont *myfont = gClient->GetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");

   cancel = new TGTextButton(cframe2, "Cancel ");
   cancel->SetText("&Cancel ");
   //if (myfont) cancel->SetFont(myfont->GetFontHandle());
   ok->SetEnabled(kTRUE);
   cancel->SetTextColor(yellow);
   cancel->SetState(kButtonEngaged);
   cframe2->AddFrame(cancel, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,3,2,2,2));

   exit = new TGTextButton(cframe2, "&Exit ","gApplication->Terminate(0)");
   cframe2->AddFrame(exit, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2,0,2,2));
   exit->SetText("&Exit ");

   AddFrame(cframe2, new TGLayoutHints(kLHintsCenterX, 2, 2, 5, 1));

   SetWindowName("Buttons' Layout");
   MapSubwindows();
   Layout();
   Resize(GetDefaultSize());
   SetWMPosition(0, 0);
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Test layout and different states of some buttons.

void testLayout()
{
   ButtonLayoutWindow *f = new ButtonLayoutWindow(gClient->GetRoot(), 100, 100);
   ProcessFrame((TGMainFrame*)f, "Buttons 1 (layout)");
   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

class TextMargin : public TGHorizontalFrame {

protected:
   TGNumberEntry *fEntry;

public:
   TextMargin(const TGWindow *p, const char *name) : TGHorizontalFrame(p)
   {
      fEntry = new TGNumberEntry(this, 0, 6, -1, TGNumberFormat::kNESInteger);
      AddFrame(fEntry, new TGLayoutHints(kLHintsLeft));
      TGLabel *label = new TGLabel(this, name);
      AddFrame(label, new TGLayoutHints(kLHintsLeft, 10));
   }
   TGTextEntry *GetEntry() const { return fEntry->GetNumberEntry(); }

};

class TextAlignWindow : public TGMainFrame {

protected:
   TGTextButton *fButton;   // button being tested

public:
   TextAlignWindow();
   void SetTextPosition(Int_t hid, Int_t vid);

};

////////////////////////////////////////////////////////////////////////////////
/// Set text position (alignment).

void TextAlignWindow::SetTextPosition(Int_t hid, Int_t vid)
{
   Int_t tj = fButton->GetTextJustify();
   tj &= ~kTextCenterX;
   tj &= ~kTextLeft;
   tj &= ~kTextRight;
   tj &= ~kTextCenterY;
   tj &= ~kTextTop;
   tj &= ~kTextBottom;
   tj |= hid;
   tj |= vid;
   fButton->SetTextJustify(tj);
}

////////////////////////////////////////////////////////////////////////////////
/// Main test window.

TextAlignWindow::TextAlignWindow() : TGMainFrame(gClient->GetRoot(), 10, 10, kHorizontalFrame)
{
   SetCleanup(kDeepCleanup);

   // Controls on right
   TGVerticalFrame *controls = new TGVerticalFrame(this);
   AddFrame(controls, new TGLayoutHints(kLHintsRight | kLHintsExpandY, 5, 5, 5, 5));

   // Separator
   TGVertical3DLine *separator = new TGVertical3DLine(this);
   AddFrame(separator, new TGLayoutHints(kLHintsRight | kLHintsExpandY));

   // Contents
   TGHorizontalFrame *contents = new TGHorizontalFrame(this);
   AddFrame(contents, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5));

   // The button for test
   fButton = new TGTextButton(contents,
                      "&This button has a multi-line label\nand shows features\navailable in the button classes");
   fButton->Resize(300, 200);
   fButton->ChangeOptions(fButton->GetOptions() | kFixedSize);
   fButton->SetToolTipText("The assigned tooltip\ncan be multi-line also", 200);
   contents->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 20, 20, 20, 20));

   TGGroupFrame *group = new TGGroupFrame(controls, "Enable/Disable");
   group->SetTitlePos(TGGroupFrame::kCenter);
   TGCheckButton *disable = new TGCheckButton(group, "Switch state\nEnable/Disable");
   disable->SetOn();
   group->AddFrame(disable, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   controls->AddFrame(group, new TGLayoutHints(kLHintsExpandX));

   // control horizontal position of the text
   TGButtonGroup *horizontal = new TGButtonGroup(controls, "Horizontal Position");
   horizontal->SetTitlePos(TGGroupFrame::kCenter);
   new TGRadioButton(horizontal, "Center", kTextCenterX);
   new TGRadioButton(horizontal, "Left", kTextLeft);
   new TGRadioButton(horizontal, "Right", kTextRight);
   horizontal->SetButton(kTextCenterX);
   controls->AddFrame(horizontal, new TGLayoutHints(kLHintsExpandX));

   // control vertical position of the text
   TGButtonGroup *vertical = new TGButtonGroup(controls, "Vertical Position");
   vertical->SetTitlePos(TGGroupFrame::kCenter);
   new TGRadioButton(vertical, "Center", kTextCenterY);
   new TGRadioButton(vertical, "Top", kTextTop);
   new TGRadioButton(vertical, "Bottom", kTextBottom);
   vertical->SetButton(kTextCenterY);
   controls->AddFrame(vertical, new TGLayoutHints(kLHintsExpandX));

   // control margins of the text
   TGGroupFrame *margins = new TGGroupFrame(controls, "Text Margins");
   margins->SetTitlePos(TGGroupFrame::kCenter);

   TextMargin *left = new TextMargin(margins, "Left");
   margins->AddFrame(left, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));

   TextMargin *right = new TextMargin(margins, "Right");
   margins->AddFrame(right, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));

   TextMargin *top = new TextMargin(margins, "Top");
   margins->AddFrame(top, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));

   TextMargin *bottom = new TextMargin(margins, "Bottom");
   margins->AddFrame(bottom, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));

   controls->AddFrame(margins, new TGLayoutHints(kLHintsExpandX));

   TGTextButton *quit = new TGTextButton(controls, "Quit");
   controls->AddFrame(quit, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 0, 0, 0, 5));

   SetWindowName("Button Test");
   MapSubwindows();
   Layout();
   Resize(GetDefaultSize());
   SetWMPosition(0, 0);
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Test different text alignments in a TGTextButton.

void testTextAlign()
{
   TextAlignWindow *f = new TextAlignWindow();
   ProcessFrame((TGMainFrame*)f, "Buttons 2 (text alignment)");
   f->SetTextPosition(kTextLeft, kTextTop);
   ProcessFrame((TGMainFrame*)f, "Buttons 3 (text alignment)");
   f->SetTextPosition(kTextCenterX, kTextTop);
   ProcessFrame((TGMainFrame*)f, "Buttons 4 (text alignment)");
   f->SetTextPosition(kTextRight, kTextTop);
   ProcessFrame((TGMainFrame*)f, "Buttons 5 (text alignment)");
   f->SetTextPosition(kTextRight, kTextCenterY);
   ProcessFrame((TGMainFrame*)f, "Buttons 6 (text alignment)");
   f->SetTextPosition(kTextRight, kTextBottom);
   ProcessFrame((TGMainFrame*)f, "Buttons 7 (text alignment)");
   f->SetTextPosition(kTextCenterX, kTextBottom);
   ProcessFrame((TGMainFrame*)f, "Buttons 8 (text alignment)");
   f->SetTextPosition(kTextLeft, kTextBottom);
   ProcessFrame((TGMainFrame*)f, "Buttons 9 (text alignment)");
   f->SetTextPosition(kTextLeft, kTextCenterY);
   ProcessFrame((TGMainFrame*)f, "Buttons 10 (text alignment)");
   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

class GroupStateWindow : public TGMainFrame {

public:
   TGTextButton        *fExit;         // Exit text button
   TGVButtonGroup      *fButtonGroup;  // Button group
   TGCheckButton       *fCheckb[4];    // Check buttons
   TGRadioButton       *fRadiob[2];    // Radio buttons

public:
   GroupStateWindow(const TGWindow *p, UInt_t w, UInt_t h);
};

////////////////////////////////////////////////////////////////////////////////
/// Main window constructor.

GroupStateWindow::GroupStateWindow(const TGWindow *p, UInt_t w, UInt_t h) :
   TGMainFrame(p, w, h)
{
   SetCleanup(kDeepCleanup);

   TGHorizontalFrame *fHL2 = new TGHorizontalFrame(this, 70, 100);
   fCheckb[0] = new TGCheckButton(fHL2, new TGHotString("Enable BG"), 100);
   fCheckb[0]->SetToolTipText("Enable/Disable the button group");
   fHL2->AddFrame(fCheckb[0], new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 1, 1, 1, 1));
   fButtonGroup = new TGVButtonGroup(fHL2, "My Button Group");
   fCheckb[1] = new TGCheckButton(fButtonGroup, new TGHotString("CB 2"), 101);
   fCheckb[2] = new TGCheckButton(fButtonGroup, new TGHotString("CB 3"), 102);
   fCheckb[3] = new TGCheckButton(fButtonGroup, new TGHotString("CB 4"), 103);
   fRadiob[0] = new TGRadioButton(fButtonGroup, new TGHotString("RB 1"), 104);
   fRadiob[1] = new TGRadioButton(fButtonGroup, new TGHotString("RB 2"), 105);
   fButtonGroup->Show();

   fHL2->AddFrame(fButtonGroup, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 1, 1, 1, 1));
   AddFrame(fHL2);

   TGHorizontalFrame *fHL3 = new TGHorizontalFrame(this, 70, 100, kFixedWidth);
   fExit = new TGTextButton(fHL3, "&Exit", 106);
   fHL3->AddFrame(fExit, new TGLayoutHints(kLHintsExpandX));
   AddFrame(fHL3, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 1, 1, 1, 1));

   //Default state
   fCheckb[0]->SetOn();
   fButtonGroup->SetState(kTRUE);

   SetWindowName("My Button Group");
   MapSubwindows();
   Layout();
   Resize(GetDefaultSize());
   SetWMPosition(0, 0);
   MapWindow();

   fButtonGroup->SetRadioButtonExclusive(kTRUE);
   fRadiob[1]->SetOn();
};

////////////////////////////////////////////////////////////////////////////////
/// Test enabled/disabled state of button group.

void testGroupState()
{
   GroupStateWindow *f = new GroupStateWindow(gClient->GetRoot(), 100, 100);
   ProcessFrame((TGMainFrame*)f, "Buttons 11 (group state)");
   f->fCheckb[0]->SetOn(kFALSE);
   f->fButtonGroup->SetState(kFALSE);
   ProcessFrame((TGMainFrame*)f, "Buttons 12 (group state)");

   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

class LabelsWindow : public TGMainFrame {

private:
   TGLabel       *fLbl1, *fLbl2, *fLbl3, *fLbl4;
public:
   LabelsWindow(const TGWindow *p, UInt_t w, UInt_t h);
   void SwitchState();
};

////////////////////////////////////////////////////////////////////////////////
/// Main window constructor.

LabelsWindow::LabelsWindow(const TGWindow *p, UInt_t w, UInt_t h) :
  TGMainFrame(p, w, h)
{
   SetCleanup(kDeepCleanup);
   // label + horizontal line
   TGGC *fTextGC;
   const TGFont *font = gClient->GetFont("-*-times-bold-r-*-*-18-*-*-*-*-*-*-*");
   if (!font)
      font = gClient->GetResourcePool()->GetDefaultFont();
   FontStruct_t labelfont = font->GetFontStruct();
   GCValues_t   gval;
   gval.fMask = kGCBackground | kGCFont | kGCForeground;
   gval.fFont = font->GetFontHandle();
   gClient->GetColorByName("yellow", gval.fBackground);
   fTextGC = gClient->GetGC(&gval, kTRUE);

   ULong_t bcolor, ycolor;
   gClient->GetColorByName("yellow", ycolor);
   gClient->GetColorByName("blue", bcolor);

   // Create a main frame
   fLbl1 = new TGLabel(this, "OwnFont & Bck/ForgrColor", fTextGC->GetGC(), labelfont, kChildFrame, bcolor);
   AddFrame(fLbl1, new TGLayoutHints(kLHintsNormal, 5, 5, 3, 4));
   fLbl1->SetTextColor(ycolor);

   fLbl2 = new TGLabel(this, "Own Font & ForegroundColor", fTextGC->GetGC(), labelfont);
   AddFrame(fLbl2,  new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   fLbl2->SetTextColor(ycolor);

   fLbl3 = new TGLabel(this, "Normal Label");
   AddFrame(fLbl3,  new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

   fLbl4 = new TGLabel(this, "Multi-line label, resized\nto 300x80 pixels",
                       fTextGC->GetGC(), labelfont, kChildFrame, bcolor);
   AddFrame(fLbl4, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   fLbl4->SetTextColor(ycolor);
   fLbl4->ChangeOptions(fLbl4->GetOptions() | kFixedSize);
   fLbl4->Resize(350, 80);

   // Create a horizontal frame containing two buttons
   TGTextButton *toggle = new TGTextButton(this, "&Toggle Labels");
   toggle->SetToolTipText("Click on the button to toggle label's state (enable/disable)");
   AddFrame(toggle, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));
   TGTextButton *exit = new TGTextButton(this, "&Exit ");
   AddFrame(exit, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));

   // Set a name to the main frame
   SetWindowName("Labels");
   MapSubwindows();
   Layout();
   Resize(GetDefaultSize());
   SetWMPosition(0, 0);
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Switch state of the labels.

void LabelsWindow::SwitchState()
{
   if (fLbl1->IsDisabled()) {
      fLbl1->Enable();
      fLbl2->Enable();
      fLbl3->Enable();
      fLbl4->Enable();
   } else {
      fLbl1->Disable();
      fLbl2->Disable();
      fLbl3->Disable();
      fLbl4->Disable();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Test different styles and the enabled/disabled state of labels.

void testLabels()
{
   LabelsWindow *f = new LabelsWindow(gClient->GetRoot(), 200, 200);
   ProcessFrame((TGMainFrame*)f, "Labels 1");
   f->SwitchState();
   ProcessFrame((TGMainFrame*)f, "Labels 2");
   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

class SplitButtonWindow : public TGMainFrame {

public:
   TGCheckButton *fCButton;
   TGCheckButton *fEButton;
   TGSplitButton *fMButton;  // Split Button
   TGPopupMenu   *fPopMenu;  // TGpopupMenu that will be attached to
                             // the button.
public:
   SplitButtonWindow(const TGWindow *p, UInt_t w, UInt_t h);

};

////////////////////////////////////////////////////////////////////////////////
/// Main window constructor.

SplitButtonWindow::SplitButtonWindow(const TGWindow *p, UInt_t w, UInt_t h) :
   TGMainFrame(p, w, h)
{
   SetCleanup(kDeepCleanup);

   TGVerticalFrame *fVL = new TGVerticalFrame(this, 100, 100);
   TGHorizontalFrame *fHL = new TGHorizontalFrame(fVL, 100, 40);

   // Create a popup menu.
   fPopMenu = new TGPopupMenu(gClient->GetRoot());
   fPopMenu->AddEntry("Button &1", 1001);
   fPopMenu->AddEntry("Button &2", 1002);
   fPopMenu->DisableEntry(1002);
   fPopMenu->AddEntry("Button &3", 1003);
   fPopMenu->AddSeparator();

   // Create a split button, the menu is adopted.
   fMButton = new TGSplitButton(fHL, new TGHotString("Button &Options"),
                                fPopMenu, 1101);

   // It is possible to add entries later
   fPopMenu->AddEntry("En&try with really really long name", 1004);
   fPopMenu->AddEntry("&Exit", 1005);

   fCButton = new TGCheckButton(fHL, new TGHotString("Split"), 1102);
   fCButton->SetState(kButtonDown);

   // Add frames to their parent for layout.
   fHL->AddFrame(fCButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY,
                                             0, 10, 0, 0));
   fEButton = new TGCheckButton(fHL, new TGHotString("Enable"), 1103);
   fEButton->SetState(kButtonDown);

   // Add frames to their parent for layout.
   fHL->AddFrame(fEButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY,
                                             0, 10, 0, 0));
   fHL->AddFrame(fMButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY));
   fVL->AddFrame(fHL, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY));
   AddFrame(fVL, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY));

   SetWindowName("SplitButton Test");
   MapSubwindows();
   Layout();
   Resize(GetDefaultSize());
   SetWMPosition(0, 0);
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Test the different configurations/states of a split button.

void testSplitButton()
{
   SplitButtonWindow *f = new SplitButtonWindow(gClient->GetRoot(),100,100);
   f->fCButton->SetState(kButtonDown);
   f->fEButton->SetState(kButtonDown);
   f->fMButton->SetState(kButtonUp);
   f->fMButton->SetSplit(kTRUE);
   ProcessFrame((TGMainFrame*)f, "Split Button 1");
   f->fCButton->SetState(kButtonUp);
   f->fEButton->SetState(kButtonDown);
   f->fMButton->SetState(kButtonUp);
   f->fMButton->SetSplit(kFALSE);
   ProcessFrame((TGMainFrame*)f, "Split Button 2");
   f->fCButton->SetState(kButtonDown);
   f->fEButton->SetState(kButtonUp);
   f->fMButton->SetState(kButtonDisabled);
   f->fMButton->SetSplit(kTRUE);
   ProcessFrame((TGMainFrame*)f, "Split Button 3");
   f->fCButton->SetState(kButtonUp);
   f->fEButton->SetState(kButtonUp);
   f->fMButton->SetState(kButtonDisabled);
   f->fMButton->SetSplit(kFALSE);
   ProcessFrame((TGMainFrame*)f, "Split Button 4");
   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

class GroupBox : public TGGroupFrame {
private:
   TGComboBox  *fCombo; // combo box
   TGTextEntry *fEntry; // text entry

public:
   GroupBox(const TGWindow *p, const char *name, const char *title);
   TGTextEntry *GetEntry() const { return fEntry; }
   TGComboBox  *GetCombo() const { return fCombo; }
};

////////////////////////////////////////////////////////////////////////////////
/// Group frame containing combobox and text entry.

GroupBox::GroupBox(const TGWindow *p, const char *name, const char *title) :
   TGGroupFrame(p, name)
{
   TGHorizontalFrame *horz = new TGHorizontalFrame(this);
   AddFrame(horz, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));
   TGLabel *label = new TGLabel(horz, title);
   horz->AddFrame(label, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   fCombo = new TGComboBox(horz);
   horz->AddFrame(fCombo, new TGLayoutHints(kLHintsRight | kLHintsExpandY,
                                            5, 0, 5, 5));
   fCombo->Resize(100, 20);

   fEntry = new TGTextEntry(this);
   AddFrame(fEntry, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));
}

////////////////////////////////////////////////////////////////////////////////

class TextEntryWindow : public TGMainFrame {

public:
   GroupBox    *fEcho;     // echo mode (echo, password, no echo)
   GroupBox    *fAlign;    // alignment (left, right, center)
   GroupBox    *fAccess;   // read-only mode
   GroupBox    *fBorder;   // border mode

public:
   TextEntryWindow();
};

////////////////////////////////////////////////////////////////////////////////
/// Main window constructor.

TextEntryWindow::TextEntryWindow() : TGMainFrame(gClient->GetRoot(), 10, 10,
                                                 kVerticalFrame)
{
   TGComboBox  *combo;
   TGTextEntry *entry;

   // recusively delete all subframes on exit
   SetCleanup(kDeepCleanup);

   fEcho = new GroupBox(this, "Echo", "Mode:");
   AddFrame(fEcho, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5, 5));
   combo = fEcho->GetCombo();
   entry = fEcho->GetEntry();
   entry->SetText("The Quick Brown Fox Jumps");
   // add entries
   combo->AddEntry("Normal", TGTextEntry::kNormal);
   combo->AddEntry("Password", TGTextEntry::kPassword);
   combo->AddEntry("No Echo", TGTextEntry::kNoEcho);
   combo->Select(TGTextEntry::kNormal);

   fAlign = new GroupBox(this, "Alignment", "Type:");
   AddFrame(fAlign, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5, 5));
   combo = fAlign->GetCombo();
   entry = fAlign->GetEntry();
   entry->SetText("Over The Lazy Dog");
   // add entries
   combo->AddEntry("Left", kTextLeft);
   combo->AddEntry("Centered", kTextCenterX);
   combo->AddEntry("Right", kTextRight);
   combo->Select(kTextLeft);

   fAccess = new GroupBox(this, "Access", "Read-only:");
   AddFrame(fAccess, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5, 5));
   combo = fAccess->GetCombo();
   entry = fAccess->GetEntry();
   entry->SetText("The Quick Brown Fox Jumps");
   // add entries
   combo->AddEntry("False", 1);
   combo->AddEntry("True", 0);
   combo->Select(1);

   fBorder = new GroupBox(this, "Border", "Drawn:");
   AddFrame(fBorder, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5, 5));
   combo = fBorder->GetCombo();
   entry = fBorder->GetEntry();
   entry->SetText("Over The Lazy Dog");
   // add entries
   combo->AddEntry("False", 0);
   combo->AddEntry("True", 1);
   combo->Select(1);

   SetWindowName("Text Entries");
   MapSubwindows();
   Layout();
   Resize(GetDefaultSize());
   SetWMPosition(0, 0);
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Test the different modes available for text entries.

void testTextEntries()
{
   TGComboBox  *combo;
   TGTextEntry *entry;

   TextEntryWindow *f = new TextEntryWindow();
   ProcessFrame((TGMainFrame*)f, "Text Entries 1");

   combo = f->fEcho->GetCombo();
   entry = f->fEcho->GetEntry();
   combo->Select(TGTextEntry::kPassword);
   entry->SetEchoMode(TGTextEntry::kPassword);

   combo = f->fAlign->GetCombo();
   entry = f->fAlign->GetEntry();
   combo->Select(kTextCenterX);
   entry->SetAlignment(kTextCenterX);

   combo = f->fAccess->GetCombo();
   entry = f->fAccess->GetEntry();
   combo->Select(0);
   entry->SetEnabled(0);

   combo = f->fBorder->GetCombo();
   entry = f->fBorder->GetEntry();
   combo->Select(0);
   entry->SetFrameDrawn(0);

   ProcessFrame((TGMainFrame*)f, "Text Entries 2");

   combo = f->fEcho->GetCombo();
   entry = f->fEcho->GetEntry();
   combo->Select(TGTextEntry::kNoEcho);
   entry->SetEchoMode(TGTextEntry::kNoEcho);

   combo = f->fAlign->GetCombo();
   entry = f->fAlign->GetEntry();
   combo->Select(kTextRight);
   entry->SetAlignment(kTextRight);

   combo = f->fAccess->GetCombo();
   entry = f->fAccess->GetEntry();
   combo->Select(1);
   entry->SetEnabled(1);

   combo = f->fBorder->GetCombo();
   entry = f->fBorder->GetEntry();
   combo->Select(1);
   entry->SetFrameDrawn(1);

   ProcessFrame((TGMainFrame*)f, "Text Entries 3");

   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

class ListTreeWindow : public TGMainFrame {

private:
   TGListTree    *fListTree;
   TGCanvas      *fCanvas;
   TGViewPort    *fViewPort;
   TList         *fNamesList;
   TGTextButton  *fTextButton;

public:
   ListTreeWindow(const TGWindow *p, UInt_t w, UInt_t h);
   void FillListTree();
   void SwitchState();
};

////////////////////////////////////////////////////////////////////////////////
/// Main window constructor.

ListTreeWindow::ListTreeWindow(const TGWindow *p, UInt_t w, UInt_t h) :
   TGMainFrame(p, w, h, kVerticalFrame)
{
   SetCleanup(kDeepCleanup);

   // canvas widget
   fCanvas = new TGCanvas(this, 300, 300);

   // canvas viewport
   fViewPort = fCanvas->GetViewPort();

   // list tree
   fListTree = new TGListTree(fCanvas, kHorizontalFrame);
   fListTree->SetCheckMode(TGListTree::kRecursive);

   fListTree->AddItem(0, "rootNode");

   fNamesList = new TList();

   TString *s1 = new TString("P1D/pclhcb08_CC-PC_HSys/CCPCAlg/Phist1");
   fNamesList->Add((TObject*)s1);

   TString *s2 = new TString("H1D/pclhcb08_CC-PC_HSys/CCPCAlg/TELL1Mult_$T3/L1/Q0/myTell");
   fNamesList->Add((TObject*)s2);

   TString *s3 = new TString("H1D/pclhcb08_CC-PC_HSys/CCPCAlg/TELL1Mult_$T3/L2/Q1/myTell");
   fNamesList->Add((TObject*)s3);

   TString *s4 = new TString("H1D/pclhcb08_CC-PC_HSys/CCPCAlg/TELL1Mult_$T3/L2/Q0/myTell");
   fNamesList->Add((TObject*)s4);

   TString *s5 = new TString("H2D/pclhcb08_CC-PC_HSys/CCPCAlg/Hist7");
   fNamesList->Add((TObject*)s5);

   TString *s6 = new TString("P1D/sdf/GaudiExample/xyProfile");
   fNamesList->Add((TObject*)s6);

   TString *s7 = new TString("H2D/sdf/GaudiExample/xyPositionPlot");
   fNamesList->Add((TObject*)s7);

   TString *s8 = new TString("H1D/sdf/GaudiExample/Mass");
   fNamesList->Add((TObject*)s8);

   TString *s9 = new TString("H1D/sdf/GaudiExample/eventtype");
   fNamesList->Add((TObject*)s9);

   AddFrame(fCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop |
                                       kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 3));
   fTextButton = new TGTextButton(this,"Text &Button...");
   AddFrame(fTextButton, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 5, 2, 5));

   FillListTree();

   SetWindowName("List Tree Test");
   MapSubwindows();
   Layout();
   Resize(490, 350);
   SetWMPosition(0, 0);
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the list tree with some hierarchical structures.

void ListTreeWindow::FillListTree()
{
   TString *listItem = 0;
   TGListTreeItem *node = 0, *histoNode = 0;

   TIter NextName(fNamesList);
   while ((listItem = (TString*)NextName())) {
      TObject *nameItem = 0;
      TObjArray *nameItems = listItem->Tokenize("/");
      TIter Next(nameItems);
      node = fListTree->GetFirstItem();
      fListTree->OpenItem(node);
      while ((nameItem = Next())) {
         if (fListTree->FindChildByName(node, nameItem->GetName())) {
            node = fListTree->FindChildByName(node, nameItem->GetName());
         } else {
            node = fListTree->AddItem(node, nameItem->GetName());
            fListTree->SetCheckBox(node, kTRUE);
            fListTree->ToggleItem(node);
            if (strcmp(nameItem->GetName(), "TELL1Mult_$T3"))
               fListTree->OpenItem(node);
            else
               fListTree->HighlightItem(node);
         }
         node->SetUserData(0);
         if (nameItem == nameItems->At(nameItems->GetEntriesFast()-2)) {
            histoNode = fListTree->AddItem(node, nameItems->Last()->GetName());
            fListTree->SetCheckBox(histoNode, kTRUE);
            fListTree->ToggleItem(histoNode);
            histoNode->SetUserData(listItem);
            break;
         }
      }
   }
   fListTree->ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Switch status of a couple of entries, to verify the propagation to
/// parent/children list tree items.

void ListTreeWindow::SwitchState()
{
   TGListTreeItem *root = 0, *node = 0;
   root = fListTree->GetFirstItem();
   node = fListTree->FindChildByName(root, "P1D");
   if (node) node = fListTree->FindChildByName(node, "pclhcb08_CC-PC_HSys");
      if (node) fListTree->CheckAllChildren(node, kTRUE);
   if (node && node->GetParent()) node->GetParent()->UpdateState();

   node = fListTree->FindChildByName(root, "H1D");
   if (node) {
      fListTree->CheckAllChildren(node, kTRUE);
      node = fListTree->FindChildByName(node, "sdf");
      if (node) node = fListTree->FindChildByName(node, "GaudiExample");
         if (node) node = fListTree->FindChildByName(node, "eventtype");
            if (node) fListTree->CheckItem(node, kFALSE);
      while (node) {
         node->UpdateState();
         node = node->GetParent();
      }
   }
   fListTree->ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Test TGListTree and TGListTreeItems.

void testListTree()
{
   ListTreeWindow *f = new ListTreeWindow(gClient->GetRoot(), 100, 100);
   ProcessFrame((TGMainFrame*)f, "List Tree 1");
   f->SwitchState();
   ProcessFrame((TGMainFrame*)f, "List Tree 2");
   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

struct shutterData_t {
   const char *pixmap_name;
   const char *tip_text;
   Int_t       id;
   TGButton   *button;
};

shutterData_t histo_data[] = {
{ "h1_s.xpm",        "TH1",      1001,  0 },
{ "h2_s.xpm",        "TH2",      1002,  0 },
{ "h3_s.xpm",        "TH3",      1003,  0 },
{ "profile_s.xpm",   "TProfile", 1004,  0 },
{ 0,                 0,          0,     0 }
};

shutterData_t function_data[] = {
{ "f1_s.xpm",        "TF1",      2001,  0 },
{ "f2_s.xpm",        "TF2",      2002,  0 },
{ 0,                 0,          0,     0 }
};

shutterData_t tree_data[] = {
{ "ntuple_s.xpm",    "TNtuple",  3001,  0 },
{ "tree_s.xpm",      "TTree",    3002,  0 },
{ "chain_s.xpm",     "TChain",   3003,  0 },
{ 0,                 0,          0,     0 }
};


class ShutterWindow : public TGMainFrame {

private:
   TGShutter        *fShutter;
   TGLayoutHints    *fLayout;
   const TGPicture  *fDefaultPic;

public:
   ShutterWindow(const TGWindow *p, UInt_t w, UInt_t h);

   void AddShutterItem(const char *name, shutterData_t *data);
   void ToggleShutterItem(const char *name);
};

////////////////////////////////////////////////////////////////////////////////
/// Create transient frame containing a shutter widget.

ShutterWindow::ShutterWindow(const TGWindow *p, UInt_t w, UInt_t h) :
   TGMainFrame(p, w, h, kVerticalFrame)
{
   // use hierarchical cleaning
   SetCleanup(kDeepCleanup);

   fDefaultPic = gClient->GetPicture("folder_s.xpm");
   fShutter = new TGShutter(this, kSunkenFrame);

   AddShutterItem("Histograms", histo_data);
   AddShutterItem("Functions", function_data);
   AddShutterItem("Trees", tree_data);

   fLayout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
   AddFrame(fShutter, fLayout);

   SetWindowName("Shutter Test");
   MapSubwindows();
   Layout();
   Resize(80, 300);
   SetWMPosition(0, 0);
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Add item in the shutter.

void ShutterWindow::AddShutterItem(const char *name, shutterData_t *data)
{
   TGShutterItem    *item;
   TGCompositeFrame *container;
   TGPictureButton  *button;
   const TGPicture  *buttonpic;
   static int id = 5001;

   TGLayoutHints *l = new TGLayoutHints(kLHintsTop | kLHintsCenterX,
                                        5, 5, 5, 0);

   item = new TGShutterItem(fShutter, new TGHotString(name), id++);
   container = (TGCompositeFrame *) item->GetContainer();

   for (int i=0; data[i].pixmap_name != 0; i++) {
      buttonpic = gClient->GetPicture(data[i].pixmap_name);
      if (!buttonpic) {
         printf("<ShutterWindow::AddShutterItem>: missing pixmap \"%s\", using default",
                data[i].pixmap_name);
         buttonpic = fDefaultPic;
      }

      button = new TGPictureButton(container, buttonpic, data[i].id);

      button->SetPicture(buttonpic);
      if (i == 0)
         button->SetState(kButtonEngaged);
      else if (i == 1)
         button->SetState(kButtonDisabled);
      else
         button->SetDisabledPicture(buttonpic);

      container->AddFrame(button, l);
      button->SetToolTipText(data[i].tip_text);
      data[i].button = button;
   }

   fShutter->AddItem(item);
}

////////////////////////////////////////////////////////////////////////////////

void ShutterWindow::ToggleShutterItem(const char *name)
{
   Long_t id = fShutter->GetItem(name)->WidgetId();
   SendMessage(fShutter, MK_MSG(kC_COMMAND, kCM_BUTTON), id, 0);
   for (int i=0; i<20;i++) {
      gSystem->ProcessEvents();
      gSystem->Sleep(10);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Test TGShutter widget.

void testShutter()
{
   ShutterWindow *f = new ShutterWindow(gClient->GetRoot(), 400, 200);
   ProcessFrame((TGMainFrame*)f, "Shutter 1");
   f->ToggleShutterItem("Functions");
   ProcessFrame((TGMainFrame*)f, "Shutter 2");
   f->ToggleShutterItem("Trees");
   ProcessFrame((TGMainFrame*)f, "Shutter 3");
   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

class ProgressbarWindow : public TGMainFrame {

private:
   TGHorizontalFrame *fHframe1;
   TGVerticalFrame   *fVframe1;
   TGLayoutHints     *fHint1, *fHint2, *fHint3, *fHint4, *fHint5;
   TGHProgressBar    *fHProg1, *fHProg2, *fHProg3;
   TGVProgressBar    *fVProg1, *fVProg2;
   TGTextButton      *fGO;

public:
   ProgressbarWindow(const TGWindow *p, UInt_t w, UInt_t h);

   void SetValues(Int_t ver);
};

////////////////////////////////////////////////////////////////////////////////
/// Main window constructor.

ProgressbarWindow::ProgressbarWindow(const TGWindow *p, UInt_t w, UInt_t h) :
   TGMainFrame(p, w, h, kHorizontalFrame)
{
   // use hierarchical cleaning
   SetCleanup(kDeepCleanup);

   fHframe1 = new TGHorizontalFrame(this, 0, 0, 0);

   fVProg1 = new TGVProgressBar(fHframe1, TGProgressBar::kFancy, 300);
   fVProg1->SetBarColor("purple");
   fVProg2 = new TGVProgressBar(fHframe1, TGProgressBar::kFancy, 300);
   fVProg2->SetFillType(TGProgressBar::kBlockFill);
   fVProg2->SetBarColor("green");

   fHframe1->Resize(300, 300);

   fVframe1 = new TGVerticalFrame(this, 0, 0, 0);

   fHProg1 = new TGHProgressBar(fVframe1, 300);
   fHProg1->ShowPosition();
   fHProg2 = new TGHProgressBar(fVframe1, TGProgressBar::kFancy, 300);
   fHProg2->SetBarColor("lightblue");
   fHProg2->ShowPosition(kTRUE, kFALSE, "%.0f events");
   fHProg3 = new TGHProgressBar(fVframe1, TGProgressBar::kStandard, 300);
   fHProg3->SetFillType(TGProgressBar::kBlockFill);

   fGO = new TGTextButton(fVframe1, "Go", 10);

   fVframe1->Resize(300, 300);

   fHint1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandY, 5, 10, 5, 5);
   fHint2 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 5, 5,  5, 10);
   fHint3 = new TGLayoutHints(kLHintsTop | kLHintsRight, 0, 50, 50, 0);
   fHint4 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandY, 0, 0, 0, 0);
   fHint5 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 0, 0);

   fHframe1->AddFrame(fVProg1, fHint1);
   fHframe1->AddFrame(fVProg2, fHint1);

   fVframe1->AddFrame(fHProg1, fHint2);
   fVframe1->AddFrame(fHProg2, fHint2);
   fVframe1->AddFrame(fHProg3, fHint2);
   fVframe1->AddFrame(fGO,     fHint3);

   AddFrame(fHframe1, fHint4);
   AddFrame(fVframe1, fHint5);

   SetWindowName("Text Entries");
   MapSubwindows();
   Layout();
   Resize(GetDefaultSize());
   SetWMPosition(0, 0);
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Set some values to our progress bars.

void ProgressbarWindow::SetValues(Int_t ver)
{
   fVProg1->Reset(); fVProg2->Reset();
   fHProg1->Reset(); fHProg2->Reset(); fHProg3->Reset();
   if (ver == 0) {
      fVProg2->SetBarColor("green");
      fVProg1->SetPosition(25.0);
      fVProg2->SetPosition(50.0);
      fHProg1->SetPosition(0.0);
      fHProg2->SetPosition(25.0);
      fHProg3->SetPosition(50.0);
   }
   else {
      fVProg2->SetBarColor("red");
      fVProg1->SetPosition(50.0);
      fVProg2->SetPosition(75.0);
      fHProg1->SetPosition(50.0);
      fHProg2->SetPosition(75.0);
      fHProg3->SetPosition(100.0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Test several styles of progress bar.

void testProgressBar()
{
   ProgressbarWindow *f = new ProgressbarWindow(gClient->GetRoot(), 600, 300);
   f->SetValues(0);
   ProcessFrame((TGMainFrame*)f, "Progress Bars 1");
   f->SetValues(1);
   ProcessFrame((TGMainFrame*)f, "Progress Bars 2");
   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

class NumberEntryWindow : public TGMainFrame {

private:
   TGVerticalFrame      *fF1;
   TGVerticalFrame      *fF2;
   TGHorizontalFrame    *fF[13];
   TGLayoutHints        *fL1;
   TGLayoutHints        *fL2;
   TGLayoutHints        *fL3;
   TGLabel              *fLabel[13];
   TGNumberEntry        *fNumericEntries[13];
   TGCheckButton        *fLowerLimit;
   TGCheckButton        *fUpperLimit;
   TGNumberEntry        *fLimits[2];
   TGCheckButton        *fPositive;
   TGCheckButton        *fNonNegative;
   TGButton             *fSetButton;
   TGButton             *fExitButton;

public:
   NumberEntryWindow(const TGWindow *p);
};

const char *numlabel[] = {
   "Integer",
   "One digit real",
   "Two digit real",
   "Three digit real",
   "Four digit real",
   "Real",
   "Degree.min.sec",
   "Min:sec",
   "Hour:min",
   "Hour:min:sec",
   "Day/month/year",
   "Month/day/year",
   "Hex"
};

const Double_t numinit[] = {
   12345, 1.0, 1.00, 1.000, 1.0000, 1.2E-12,
   90 * 3600, 120 * 60, 12 * 60, 12 * 3600 + 15 * 60,
   19991121, 19991121, (Double_t) 0xDEADFACEU
};

////////////////////////////////////////////////////////////////////////////////
/// Main window constructor.

NumberEntryWindow::NumberEntryWindow(const TGWindow *p) :
   TGMainFrame(p, 10, 10, kHorizontalFrame)
{
   // use hierarchical cleaning
   SetCleanup(kDeepCleanup);

   TGGC myGC = *gClient->GetResourcePool()->GetFrameGC();
   TGFont *myfont = gClient->GetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
   if (myfont) myGC.SetFont(myfont->GetFontHandle());

   fF1 = new TGVerticalFrame(this, 200, 300);
   fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2);
   AddFrame(fF1, fL1);
   fL2 = new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2);
   for (int i = 0; i < 13; i++) {
      fF[i] = new TGHorizontalFrame(fF1, 200, 30);
      fF1->AddFrame(fF[i], fL2);
      fNumericEntries[i] = new TGNumberEntry(fF[i], numinit[i], 12, i + 20,
                                             (TGNumberFormat::EStyle) i);
      fF[i]->AddFrame(fNumericEntries[i], fL2);
      fLabel[i] = new TGLabel(fF[i], numlabel[i], myGC(), myfont->GetFontStruct());
      fF[i]->AddFrame(fLabel[i], fL2);
   }
   fF2 = new TGVerticalFrame(this, 200, 500);
   fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2);
   AddFrame(fF2, fL3);
   fLowerLimit = new TGCheckButton(fF2, "lower limit:", 4);
   fF2->AddFrame(fLowerLimit, fL3);
   fLimits[0] = new TGNumberEntry(fF2, 0, 12, 10);
   fLimits[0]->SetLogStep(kFALSE);
   fF2->AddFrame(fLimits[0], fL3);
   fUpperLimit = new TGCheckButton(fF2, "upper limit:", 5);
   fF2->AddFrame(fUpperLimit, fL3);
   fLimits[1] = new TGNumberEntry(fF2, 0, 12, 11);
   fLimits[1]->SetLogStep(kFALSE);
   fF2->AddFrame(fLimits[1], fL3);
   fPositive = new TGCheckButton(fF2, "Positive", 6);
   fF2->AddFrame(fPositive, fL3);
   fNonNegative = new TGCheckButton(fF2, "Non negative", 7);
   fF2->AddFrame(fNonNegative, fL3);
   fSetButton = new TGTextButton(fF2, " Set ", 2);
   fF2->AddFrame(fSetButton, fL3);
   fExitButton = new TGTextButton(fF2, " Close ", 1);
   fF2->AddFrame(fExitButton, fL3);

   // set dialog box title
   SetWindowName("Number Entry Test");
   MapSubwindows();
   Layout();
   Resize(GetDefaultSize());
   SetWMPosition(0, 0);
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Test number entries in different formats.

void testNumberEntry()
{
   NumberEntryWindow *f = new NumberEntryWindow(gClient->GetRoot());
   ProcessFrame((TGMainFrame*)f, "Number Entries");
   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

const char *editortxt1 =
"This is the ROOT text edit widget TGTextEdit. It is not intended as\n"
"a full developers editor, but it is relatively complete and can ideally\n"
"be used to edit scripts or to present users editable config files, etc.\n\n"
"The text edit widget supports standard emacs style ctrl-key navigation\n"
"in addition to the arrow keys. By default the widget has under the right\n"
"mouse button a popup menu giving access to several built-in functions.\n\n"
"Cut, copy and paste between different editor windows and any other\n"
"standard text handling application is supported.\n\n"
"Text can be selected with the mouse while holding the left button\n"
"or with the arrow keys while holding the shift key pressed. Use the\n"
"middle mouse button to paste text at the current mouse location."
;
const char *editortxt2 =
"Mice with scroll-ball are properly supported.\n\n"
"This are the currently defined key bindings:\n"
"Left Arrow\n"
"    Move the cursor one character leftwards.\n"
"    Scroll when cursor is out of frame.\n"
"Right Arrow\n"
"    Move the cursor one character rightwards.\n"
"    Scroll when cursor is out of frame.\n"
"Backspace\n"
"    Deletes the character on the left side of the text cursor and moves the\n"
"    cursor one position to the left. If a text has been marked by the user"
;
const char *editortxt3 =
"    (e.g. by clicking and dragging) the cursor will be put at the beginning\n"
"    of the marked text and the marked text will be removed.\n"
"Home\n"
"    Moves the text cursor to the left end of the line. If mark is TRUE text\n"
"    will be marked towards the first position, if not any marked text will\n"
"    be unmarked if the cursor is moved.\n"
"End\n"
"    Moves the text cursor to the right end of the line. If mark is TRUE text\n"
"    will be marked towards the last position, if not any marked text will\n"
"    be unmarked if the cursor is moved.\n"
;

class EditorWindow : public TGMainFrame {

private:
   TGTextEdit       *fEdit;   // text edit widget
   TGTextButton     *fOK;     // OK button
   TGLayoutHints    *fL1;     // layout of TGTextEdit
   TGLayoutHints    *fL2;     // layout of OK button

public:
   EditorWindow(const TGWindow *p, UInt_t w, UInt_t h);

   void   LoadBuffer(const char *buffer);
   void   AddBuffer(const char *buffer);

};

////////////////////////////////////////////////////////////////////////////////
/// Main window constructor.

EditorWindow::EditorWindow(const TGWindow *p, UInt_t w, UInt_t h) :
   TGMainFrame(p, w, h)
{
   // use hierarchical cleaning
   SetCleanup(kDeepCleanup);

   fEdit = new TGTextEdit(this, w, h, kSunkenFrame | kDoubleBorder);
   fL1 = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 3, 3, 3, 3);
   AddFrame(fEdit, fL1);

   // set selected text colors
   Pixel_t pxl;
   gClient->GetColorByName("#ccccff", pxl);
   fEdit->SetSelectBack(pxl);
   fEdit->SetSelectFore(TGFrame::GetBlackPixel());

   fOK = new TGTextButton(this, "  &OK  ");
   fL2 = new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 0, 0, 5, 5);
   AddFrame(fOK, fL2);

   SetWindowName("Number Entry Test");
   MapSubwindows();
   Layout();
   Resize(GetDefaultSize());
   SetWMPosition(0, 0);
   LoadBuffer(editortxt1);
   AddBuffer(editortxt2);
   AddBuffer(editortxt3);
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Load a text buffer in the editor.

void EditorWindow::LoadBuffer(const char *buffer)
{
   fEdit->LoadBuffer(buffer);
}

////////////////////////////////////////////////////////////////////////////////
/// Add text to the editor.

void EditorWindow::AddBuffer(const  char *buffer)
{
   TGText txt;
   txt.LoadBuffer(buffer);
   fEdit->AddText(&txt);
}

////////////////////////////////////////////////////////////////////////////////
/// Very simple test of the TGTextEdit widget.

void testEditor()
{
   EditorWindow *f = new EditorWindow(gClient->GetRoot(), 550, 400);
   ProcessFrame((TGMainFrame*)f, "Text Editor");
   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Test Simple Canvas...

void testCanvas()
{
   TCanvas *c = new TCanvas("c", "Test Canvas", 0, 0, 800, 600);
   TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
   ProcessFrame((TGMainFrame*)rc, "ROOT Canvas");
   c->Close();
   delete c;
}

////////////////////////////////////////////////////////////////////////////////

class MyColorDialog : public TGColorDialog {

public:
   MyColorDialog(const TGWindow *p = 0, const TGWindow *m = 0, Int_t *retc = 0,
                 Pixel_t *color = 0);

   void   Close() { CloseWindow(); }
   void   SwitchTab(Int_t pos);
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

MyColorDialog::MyColorDialog(const TGWindow *p, const TGWindow *m, Int_t *retc,
                             ULong_t *color) : TGColorDialog(p, m, retc, color,
                             kFALSE)
{
   SetWMPosition(0, 0);
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////

void MyColorDialog::SwitchTab(Int_t pos)
{
   fTab->SetTab(pos, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Test Color Selection Dialog.

void testColorDlg()
{
   Int_t retc = 0;
   ULong_t color = 0xcfcfcf;
   MyColorDialog *dlg = new MyColorDialog(gClient->GetDefaultRoot(),
                                          gClient->GetDefaultRoot(),
                                          &retc, &color);
   ProcessFrame((TGMainFrame*)dlg, "Color Dialog 1");
   dlg->SwitchTab(1);
   ProcessFrame((TGMainFrame*)dlg, "Color Dialog 2");
   dlg->Close();
   delete dlg;
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Test Font Selection Dialog.

void testFontDlg()
{
   TGFontDialog::FontProp_t prop;
   TGFontDialog *dlg = new TGFontDialog(gClient->GetDefaultRoot(),
                                        gClient->GetDefaultRoot(),
                                        &prop, "", 0, kFALSE);
   dlg->SetWMPosition(0, 0);
   dlg->MapWindow();
   gSystem->Sleep(250);
   ProcessFrame((TGMainFrame*)dlg, "Font Dialog");
   delete dlg;
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Test Search Dialog.

void testSearchDlg()
{
   Int_t retc = 0;
   TGSearchType sstruct;
   sstruct.fClose = kFALSE;
   TGSearchDialog *dlg = new TGSearchDialog(gClient->GetDefaultRoot(),
                                            gClient->GetDefaultRoot(),
                                            1, 1, &sstruct, &retc);
   dlg->SetWMPosition(0, 0);
   dlg->MapWindow();
   ProcessFrame((TGMainFrame*)dlg, "Search Dialog");
   delete dlg;
}

////////////////////////////////////////////////////////////////////////////////

class LabelEntry : public TGCompositeFrame {
public:
   TGLabel * fLabel;
   TGNumberEntry * fNumberEntry;
   TGTextButton * fTextButton1, * fTextButton2;
   TGTableLayout * tabLayout;

public:
   LabelEntry(const TGWindow * parent, Int_t Width, Int_t Height);
};

////////////////////////////////////////////////////////////////////////////////

LabelEntry::LabelEntry(const TGWindow * parent, Int_t Width, Int_t Height) :
   TGCompositeFrame(parent, Width, Height)
{
   tabLayout=new TGTableLayout(this, 2,5, kFALSE);
   SetLayoutManager(tabLayout);

   fLabel = new TGLabel(this, "Label");
   fLabel->ChangeOptions(fLabel->GetOptions() | kFixedWidth);
   fNumberEntry = new TGNumberEntry(this, 0., 30, 1);
   fTextButton1 = new TGTextButton(this, "TextButton");
   fTextButton2 = new TGTextButton(this, "TextButton2");

   AddFrame(fLabel, new TGTableLayoutHints(0,1,0,2,
            kLHintsExpandX | kLHintsExpandY |
            kLHintsCenterX | kLHintsCenterY |
            kLHintsFillX | kLHintsFillY,
            1, 1, 1, 1));
   AddFrame(fNumberEntry, new TGTableLayoutHints(1,2,0,2,
            kLHintsExpandX | kLHintsExpandY |
            kLHintsCenterX | kLHintsCenterY |
            kLHintsFillX|kLHintsFillY,
            1, 1, 1, 1));
   AddFrame(fTextButton1, new TGTableLayoutHints(2,3,0,1,
            kLHintsExpandX | kLHintsExpandY |
            kLHintsCenterX | kLHintsCenterY |
            kLHintsFillX | kLHintsFillY,
            1, 1, 1, 1));
   AddFrame(fTextButton2, new TGTableLayoutHints(2,3,1,2,
            kLHintsExpandX | kLHintsExpandY |
            kLHintsCenterX | kLHintsCenterY |
            kLHintsFillX | kLHintsFillY,
            1, 1, 1, 1));

   Layout();
   Resize(Width, Height);
   ChangeOptions(GetOptions() | kFixedWidth);
}

////////////////////////////////////////////////////////////////////////////////
/// Test table layout.

void testTableLayout()
{
   LabelEntry *l1, *l2;
   TGMainFrame *mf = new TGMainFrame (gClient->GetRoot(), 700, 200);
   mf->SetCleanup(kDeepCleanup);

   l1= new LabelEntry(mf, 500, 60);
   l2= new LabelEntry(mf, 600, 70);
   mf->AddFrame(l1);
   mf->AddFrame(l2);
   mf->SetWindowName("Test Table Layout");
   mf->MapSubwindows();
   mf->Layout();
   mf->Resize(600, 130);
   mf->SetWMPosition(0, 0);
   mf->MapWindow();
   ProcessFrame((TGMainFrame*)mf, "Table Layout");
   mf->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Test the TGPack widget.

void testPack()
{
   TGPack *hp = 0;
   TGPack *vp = 0;
   TGTextButton* b = 0;

   TGMainFrame* mf = new TGMainFrame(0, 400, 300);
   mf->SetCleanup(kDeepCleanup);
   mf->SetWindowName("Foo");

   hp = new TGPack(mf, mf->GetWidth(), mf->GetHeight());
   hp->SetVertical(kFALSE);

   b = new TGTextButton(hp, "Ailaaha");  hp->AddFrame(b);

   vp = new TGPack(hp, hp->GetWidth(), hp->GetHeight());
   b = new TGTextButton(vp, "Blaaaaa");  vp->AddFrameWithWeight(b, 0, 5);
   b = new TGTextButton(vp, "Blooooo");  vp->AddFrameWithWeight(b, 0, 3);
   b = new TGTextButton(vp, "Bleeeee");  vp->AddFrameWithWeight(b, 0, 5);
   hp->AddFrame(vp, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   b = new TGTextButton(hp, "Cilnouk");  hp->AddFrame(b);

   mf->AddFrame(hp, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   mf->MapSubwindows();
   mf->Layout();
   mf->SetWMPosition(0, 0);
   mf->MapWindow();
   gSystem->ProcessEvents();
   gSystem->Sleep(10);
   ProcessFrame((TGMainFrame*)mf, "Pack Frames");
   mf->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

enum ESlidersIds {
   VSId1,
   HSId1,
   VSId2,
   HSId2,
   VSId3,
   HSId3
};

class SliderWindow : public TGMainFrame {

private:
   TGVerticalFrame   *fVframe1, *fVframe2, *fVframe3;
   TGLayoutHints     *fBly, *fBfly1;
   TGHSlider         *fHslider1;
   TGVSlider         *fVslider1;
   TGDoubleVSlider   *fVslider2;
   TGDoubleHSlider   *fHslider2;
   TGTripleVSlider   *fVslider3;
   TGTripleHSlider   *fHslider3;

public:
   SliderWindow(const TGWindow *p, UInt_t w, UInt_t h);
};

SliderWindow::SliderWindow(const TGWindow *p, UInt_t w, UInt_t h) :
   TGMainFrame(p, w, h)
{
   // Dialog used to test the different supported sliders.

   // use hierarchical cleaning
   SetCleanup(kDeepCleanup);

   ChangeOptions((GetOptions() & ~kVerticalFrame) | kHorizontalFrame);

   fVframe1 = new TGVerticalFrame(this, 0, 0, 0);

   fHslider1 = new TGHSlider(fVframe1, 100, kSlider1 | kScaleBoth, HSId1);
   fHslider1->SetRange(0,50);

   fVslider1 = new TGVSlider(fVframe1, 100, kSlider2 | kScaleBoth, VSId1);
   fVslider1->SetRange(0,8);

   fVframe1->Resize(100, 100);

   fVframe2 = new TGVerticalFrame(this, 0, 0, 0);

   fHslider2 = new TGDoubleHSlider(fVframe2, 150, kDoubleScaleBoth, HSId2);
   fHslider2->SetRange(0,3);

   fVslider2 = new TGDoubleVSlider(fVframe2, 100, kDoubleScaleBoth, VSId2);

   fVslider2->SetRange(-10,10);
   fVframe2->Resize(100, 100);

   fVframe3 = new TGVerticalFrame(this, 0, 0, 0);

   fHslider3 = new TGTripleHSlider(fVframe3, 100, kDoubleScaleBoth, HSId3,
                                   kHorizontalFrame, GetDefaultFrameBackground(),
                                   kFALSE, kFALSE, kFALSE, kFALSE);
   fHslider3->SetRange(0.05,5.0);

   fVslider3 = new TGTripleVSlider(fVframe3, 100, kDoubleScaleBoth, VSId3,
                                   kVerticalFrame, GetDefaultFrameBackground(),
                                   kFALSE, kFALSE, kFALSE, kFALSE);
   fVslider3->SetRange(0.05,5.0);
   fVframe3->Resize(100, 100);

   //--- layout for buttons: top align, equally expand horizontally
   fBly = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 3, 0);

   //--- layout for the frame: place at bottom, right aligned
   fBfly1 = new TGLayoutHints(kLHintsTop | kLHintsRight, 20, 10, 15, 0);

   fVframe1->AddFrame(fHslider1, fBly);
   fVframe1->AddFrame(fVslider1, fBly);

   fVframe2->AddFrame(fHslider2, fBly);
   fVframe2->AddFrame(fVslider2, fBly);

   fVframe3->AddFrame(fHslider3, fBly);
   fVframe3->AddFrame(fVslider3, fBly);

   AddFrame(fVframe3, fBfly1);
   AddFrame(fVframe2, fBfly1);
   AddFrame(fVframe1, fBfly1);

   SetWindowName("Slider Test");
   TGDimension size = GetDefaultSize();
   Resize(size);

   fHslider3->SetPosition(0.15,1.5);
   fHslider3->SetPointerPosition(0.75);

   fVslider3->SetPosition(0.05,2.5);
   fVslider3->SetPointerPosition(1.0);

   SetWMSize(size.fWidth, size.fHeight);
   SetWMSizeHints(size.fWidth, size.fHeight, size.fWidth, size.fHeight, 0, 0);
   SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
                                     kMWMDecorMinimize | kMWMDecorMenu,
                      kMWMFuncAll |  kMWMFuncResize    | kMWMFuncMaximize |
                                     kMWMFuncMinimize,
                      kMWMInputModeless);

   MapSubwindows();
   MapWindow();

}

////////////////////////////////////////////////////////////////////////////////
/// Test horizontal and vertical sliders.

void testSliders()
{
   SliderWindow *f = new SliderWindow(gClient->GetRoot(), 400, 200);
   ProcessFrame((TGMainFrame*)f, "Sliders 1");
   f->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Popup the GUI...

void testBrowsers()
{
   gSystem->RedirectOutput(gTmpfilename.Data(), "w", &gRH);
   gEnv->SetValue("Browser.Name", "TRootBrowserLite");
   TBrowser *b = new TBrowser();
   gSystem->RedirectOutput(0, 0, &gRH);
   TGMainFrame *f = (TGMainFrame *)b->GetBrowserImp()->GetMainFrame();
   if (f) {
      ProcessFrame((TGMainFrame*)f, "Root Browser Lite");
      f->CloseWindow();
   }
   gSystem->RedirectOutput(gTmpfilename.Data(), "w", &gRH);
   gEnv->SetValue("Browser.Name", "TRootBrowser");
   b = new TBrowser();
   gSystem->RedirectOutput(0, 0, &gRH);
   f = (TGMainFrame *)b->GetBrowserImp()->GetMainFrame();
   if (f) {
      ProcessFrame((TGMainFrame*)f, "Root Browser");
      f->CloseWindow();
   }
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Test TGSplitFrame.

void testSplitFrame()
{
   TGMainFrame *mf = new TGMainFrame(gClient->GetDefaultRoot(), 200, 200);
   mf->SetCleanup(kDeepCleanup);
   TGSplitFrame *first = new TGSplitFrame(mf, 200, 200);
   mf->AddFrame(first, new TGLayoutHints(kLHintsExpandX | kLHintsExpandX, 0, 0, 0, 0));
   first->HSplit();
   first->GetFirst()->VSplit();
   first->GetSecond()->VSplit();
   first->GetSecond()->GetSecond()->SetEditable();
   new TGTextEditor("stressGUI.cxx", gClient->GetRoot());
   first->GetSecond()->GetSecond()->SetEditable(kFALSE);
   mf->MapSubwindows();
   mf->Resize(600, 400);
   mf->SetWMPosition(0, 0);
   mf->MapWindow();
   ProcessFrame(mf, "Split Frame 1");
   mf->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Test the ROOT control bar.

void testControlBars()
{
   TControlBar *bar = new TControlBar("vertical", "Demos",10,10);
   bar->AddButton("Help Demos",".x demoshelp.C",        "Click Here For Help on Running the Demos");
   bar->AddButton("browser",   "new TBrowser;",         "Start the ROOT Browser");
   bar->AddButton("framework", ".x graphics/framework.C","An Example of Object Oriented User Interface");
   bar->AddButton("first",     ".x graphics/first.C",   "An Example of Slide with Root");
   bar->AddButton("hsimple",   ".x hsimple.C",          "An Example Creating Histograms/Ntuples on File");
   bar->AddButton("hsum",      ".x hist/hsum.C",        "Filling Histograms and Some Graphics Options");
   bar->AddButton("formula1",  ".x graphics/formula1.C","Simple Formula and Functions");
   bar->AddButton("surfaces",  ".x graphs/surfaces.C",  "Surface Drawing Options");
   bar->AddButton("fillrandom",".x hist/fillrandom.C",  "Histograms with Random Numbers from a Function");
   bar->AddButton("fit1",      ".x fit/fit1.C",         "A Simple Fitting Example");
   bar->AddButton("multifit",  ".x fit/multifit.C",     "Fitting in Subranges of Histograms");
   bar->AddButton("h1draw",    ".x hist/h1draw.C",      "Drawing Options for 1D Histograms");
   bar->AddButton("graph",     ".x graphs/graph.C",     "Example of a Simple Graph");
   bar->AddButton("gerrors",   ".x graphs/gerrors.C",   "Example of a Graph with Error Bars");
   bar->AddButton("tornado",   ".x graphics/tornado.C", "Examples of 3-D PolyMarkers");
   bar->AddButton("shapes",    ".x geom/shapes.C",      "The Geometry Shapes");
   bar->AddButton("geometry",  ".x geom/geometry.C",    "Creation of the NA49 Geometry File");
   bar->AddButton("na49view",  ".x geom/na49view.C",    "Two Views of the NA49 Detector Geometry");
   bar->AddButton("file",      ".x io/file.C",          "The ROOT File Format");
   bar->AddButton("fildir",    ".x io/fildir.C",        "The ROOT File, Directories and Keys");
   bar->AddButton("tree",      ".x tree/tree.C",        "The Tree Data Structure");
   bar->AddButton("ntuple1",   ".x tree/ntuple1.C",     "Ntuples and Selections");
   bar->AddButton("rootmarks", ".x rootmarks.C",        "Prints an Estimated ROOTMARKS for Your Machine");
   bar->SetButtonWidth(90);
   bar->Show();
   TControlBarImp *imp = bar->GetControlBarImp();
   TGMainFrame *f = dynamic_cast<TGMainFrame*>(imp);
   if (f) {
      ProcessFrame(f, "Control Bar");
      f->CloseWindow();
   }
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Very simple test of the ROOT help dialog.

void testHelpDialog()
{
   TRootHelpDialog *hd = new TRootHelpDialog(gClient->GetRoot(), "About ROOT...", 600, 400);
   hd->SetText(gHelpAbout);
   hd->Popup();
   ProcessFrame((TGMainFrame *)hd, "Help Dialog");
   hd->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Test the ASImage palette editor.

void testPaletteEditor()
{
   const char *fname = "galaxy.root";
   TFile *gal = 0;
   if (!gSystem->AccessPathName(fname)) {
      gal = TFile::Open(fname);
   } else {
      gal = TFile::Open(Form("http://root.cern.ch/files/%s",fname));
   }
   if (!gal) return;
   TImage *img = (TImage*)gal->Get("n4254");
   //img->Draw();

   TASPaletteEditor *f = new TASPaletteEditor((TAttImage *)img, 80, 25);
   ProcessFrame((TGMainFrame*)f, "Palette Editor");
   f->CloseWindow();

   delete img;
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Test the HTML Browser.

void testHtmlBrowser()
{
   TGHtmlBrowser *b = new TGHtmlBrowser("http://bellenot.web.cern.ch/bellenot/Public/html_test/html_test.html");
   ProcessFrame((TGMainFrame*)b, "HTML Browser 1");
   b->Selected("http://bellenot.web.cern.ch/bellenot/Public/html_test/gallery/");
   ProcessFrame((TGMainFrame*)b, "HTML Browser 2");
   b->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////
//                             ROOT GUI tutorials
////////////////////////////////////////////////////////////////////////////////

// list of excluded macros
const char *excluded[] = {
   "_playback",
   "calendar",
   "customContextMenu.C",
   "customTH1Fmenu",
   "exec_macro",
   "guitest0",
   0
};

////////////////////////////////////////////////////////////////////////////////
/// start a new ROOT process for the execution of the macro.

Int_t bexec(TString &dir, const char *macro)
{
#ifdef WIN32
   return gSystem->Exec(TString::Format("set ROOT_HIST=0 & root.exe -l -q exec_macro.C(\\\"%s/%s\\\") >nul 2>&1", dir.Data(), macro));
#else
   return gSystem->Exec(TString::Format("ROOT_HIST=0 root.exe -l -q exec_macro.C\\(\\\"%s/%s\\\"\\) >&/dev/null", dir.Data(), macro));
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Run the macros available in $ROOTSYS/tutorials/gui

void run_tutorials()
{
   gClient->HandleInput();
   gSystem->Sleep(50);
   gSystem->ProcessEvents();
   TString dir = gRootSys + "/tutorials/gui";
   TString savdir = gSystem->WorkingDirectory();
   TSystemDirectory sysdir(dir.Data(), dir.Data());
   TList *files = sysdir.GetListOfFiles();

   gSystem->ChangeDirectory(dir.Data());
   dir = gRootSys + "/tutorials";
   TString reqfile = dir + "/hsimple.root";
   if (gSystem->AccessPathName(reqfile, kFileExists)) {
      bexec(dir, "hsimple.C");
      gSystem->Unlink("hsimple_1.png");
   }
   dir += "/tree";
   reqfile = dir + "/cernstaff.root";
   if (gSystem->AccessPathName(reqfile, kFileExists)) {
      bexec(dir, "cernbuild.C");
   }
   dir = gRootSys + "/tutorials/gui";

   if (files) {
      TIter next(files);
      TSystemFile *file;
      TString fname;

      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.EndsWith(".C")) {
            Bool_t skip = kFALSE;
            for (int i=0; excluded[i]; i++) {
               if (strstr(fname, excluded[i])) {
                  skip = kTRUE;
                  break;
               }
            }
            if (!skip) {
               bexec(dir, fname.Data());
               ProcessMacro(fname.Data(), fname.Data());
            }
         }
      }
      delete files;
   }
   dir = gRootSys + "/tutorials/image";
   bexec(dir, "galaxy_image.C");
   ProcessMacro("galaxy_image.C", "galaxy_image.C");

   gSystem->ChangeDirectory(savdir.Data());
}

////////////////////////////////////////////////////////////////////////////////
//                             Recorder sessions
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// count characters in the file, skipping cr/lf

Int_t file_size(const char *filename)
{
   FILE *lunin;
   Int_t c, wc = 0;

   lunin = fopen(filename, "rb");
   if (lunin == 0) return -1;
   while (!feof(lunin)) {
      c = fgetc(lunin);
      if (c != 0x0d && c != 0x0a)
         wc++;
   }
   fclose(lunin);
   return wc;
}

////////////////////////////////////////////////////////////////////////////////

void guitest_playback()
{
   Int_t i;
   Bool_t ret;
   TRecorder r;
   Int_t guitest_ref[11], guitest_err[11], guitest_size[11];

   printf("Guitest Playback..............................................");
   gSystem->RedirectOutput(gTmpfilename.Data(), "w", &gRH);
   TString savdir = gSystem->WorkingDirectory();
   TString dir = gRootSys + "/tutorials/gui";
   gSystem->ChangeDirectory(dir.Data());

   // first delete old files, if any
   for (i=0;i<11;++i) {
      gSystem->Unlink(TString::Format("%s/guitest%03d.C", dir.Data(), i+1));
   }
   TStopwatch sw;
   ret = r.Replay("http://root.cern.ch/files/guitest_playback.root");

   // wait for the recorder to finish the replay
   while (ret && r.GetState() == TRecorder::kReplaying) {
      gSystem->ProcessEvents();
      gSystem->Sleep(1);
      if (sw.RealTime() > 600.0) {
         r.ReplayStop();
         ret = kFALSE;
      }
      sw.Continue();
   }

   gSystem->RedirectOutput(0, 0, &gRH);
   if (ret) printf("... DONE\n");
   else printf(". FAILED\n");
   CloseMainframes();
   for (i=0;i<11;++i) {
      guitest_ref[i] = 0;
      guitest_err[i] = 100;
      guitest_size[i] = file_size(TString::Format("%s/guitest%03d.C",
                                  dir.Data(), i+1));
   }
   guitest_ref[0]  = 24957;
   guitest_ref[1]  = 6913;
   guitest_ref[2]  = 16402;
   guitest_ref[3]  = 10688;
   guitest_ref[4]  = 6554;
   guitest_ref[5]  = 24239;
   guitest_ref[6]  = 25069;
   guitest_ref[7]  = 25126;
   guitest_ref[8]  = 25175;
   guitest_ref[9]  = 25324;
   guitest_ref[10] = 68657;
   for (i=0;i<11;++i) {
      printf("guitest %02d: output............................................", i+1);
      if (TMath::Abs(guitest_ref[i] - guitest_size[i]) <= guitest_err[i]) {
         printf("..... OK\n");
         // delete successful tests, keep only the failing ones (for verification)
         gSystem->Unlink(TString::Format("%s/guitest%03d.C", dir.Data(), i+1));
      }
      else {
         printf(". FAILED\n");
         printf("         File Size = %d\n", guitest_size[i]);
         printf("          Ref Size = %d\n", guitest_ref[i]);
      }
   }
   gSystem->ChangeDirectory(savdir.Data());
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

void dnd_playback()
{
   Bool_t ret;
   TRecorder r;
   printf("Drag and Drop Playback........................................");
   gSystem->RedirectOutput(gTmpfilename.Data(), "w", &gRH);
   TString savdir = gSystem->WorkingDirectory();
   TString filename = "http://bellenot.web.cern.ch/bellenot/recorder/dnd_playback";
#ifdef WIN32
   filename += "_win.root";
#else
   filename += "_x11.root";
#endif
   TString dir = gRootSys + "/tutorials/gui";
   gSystem->ChangeDirectory(dir.Data());

   TStopwatch sw;
   ret = r.Replay(filename.Data());

   // wait for the recorder to finish the replay
   while (ret && r.GetState() == TRecorder::kReplaying) {
      gSystem->ProcessEvents();
      gSystem->Sleep(1);
      // add timeout...
      if (sw.RealTime() > 180.0) {
         r.ReplayStop();
         ret = kFALSE;
      }
      sw.Continue();
   }
   gSystem->RedirectOutput(0, 0, &gRH);
   if (ret) printf("... DONE\n");
   else printf(". FAILED\n");
   CloseMainframes();
   gSystem->ChangeDirectory(savdir.Data());
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

void mditest_playback()
{
   Bool_t ret;
   TRecorder r;
   printf("MDI Test Playback.............................................");
   gSystem->RedirectOutput(gTmpfilename.Data(), "w", &gRH);
   TString savdir = gSystem->WorkingDirectory();
   TString dir = gRootSys + "/tutorials/gui";
   gSystem->ChangeDirectory(dir.Data());

   TStopwatch sw;
   ret = r.Replay("http://bellenot.web.cern.ch/bellenot/recorder/mditest_playback.root");

   // wait for the recorder to finish the replay
   while (ret && r.GetState() == TRecorder::kReplaying) {
      gSystem->ProcessEvents();
      gSystem->Sleep(1);
      // add timeout...
      if (sw.RealTime() > 180.0) {
         r.ReplayStop();
         ret = kFALSE;
      }
      sw.Continue();
   }
   gSystem->RedirectOutput(0, 0, &gRH);
   if (ret) printf("... DONE\n");
   else printf(". FAILED\n");
   CloseMainframes();
   gSystem->ChangeDirectory(savdir.Data());
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

void graph_edit_playback()
{
   Bool_t ret;
   TRecorder r;
   printf("Graphic Editors Playback......................................");
   gSystem->RedirectOutput(gTmpfilename.Data(), "w", &gRH);
   TString savdir = gSystem->WorkingDirectory();
   TString dir = gRootSys + "/tutorials/graphics";
   gSystem->ChangeDirectory(dir.Data());

   TStopwatch sw;
   ret = r.Replay("http://root.cern.ch/files/graphedit_playback.root");

   // wait for the recorder to finish the replay
   while (ret && r.GetState() == TRecorder::kReplaying) {
      gSystem->ProcessEvents();
      gSystem->Sleep(1);
      // add timeout...
      if (sw.RealTime() > 300.0) {
         r.ReplayStop();
         ret = kFALSE;
      }
      sw.Continue();
   }
   gSystem->RedirectOutput(0, 0, &gRH);
   if (ret) printf("... DONE\n");
   else printf(". FAILED\n");
   CloseMainframes();
   gSystem->ChangeDirectory(savdir.Data());
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

void fitpanel_playback()
{
   Bool_t ret;
   TRecorder r;
   printf("Fit Panel Playback ...........................................");
   gSystem->RedirectOutput(gTmpfilename.Data(), "w", &gRH);
   TString savdir = gSystem->WorkingDirectory();
   TString dir = gRootSys + "/tutorials/fit";
   gSystem->ChangeDirectory(dir.Data());

   TStopwatch sw;
   ret = r.Replay("http://root.cern.ch/files/fitpanel_playback.root");

   // wait for the recorder to finish the replay
   while (ret && r.GetState() == TRecorder::kReplaying) {
      gSystem->ProcessEvents();
      gSystem->Sleep(1);
      // add timeout...
      if (sw.RealTime() > 300.0) {
         r.ReplayStop();
         ret = kFALSE;
      }
      sw.Continue();
   }
   gSystem->RedirectOutput(0, 0, &gRH);
   if (ret) printf("... DONE\n");
   else printf(". FAILED\n");
   CloseMainframes();
   gSystem->ChangeDirectory(savdir.Data());
}

