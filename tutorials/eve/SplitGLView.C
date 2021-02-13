/// \file
/// \ingroup tutorial_eve
/// Helper classes for the alice_esd_split.C demo.
///
/// \macro_code
///
/// \author Bertrand Bellenot

#include "TApplication.h"
#include "TSystem.h"
#include "TGFrame.h"
#include "TGLayout.h"
#include "TGSplitter.h"
#include "TGLWidget.h"
#include "TEvePad.h"
#include "TGeoManager.h"
#include "TString.h"
#include "TGMenu.h"
#include "TGStatusBar.h"
#include "TGFileDialog.h"
#include "TGMsgBox.h"
#include "TGLPhysicalShape.h"
#include "TGLLogicalShape.h"
#include "HelpText.h"
#include "TClass.h"
#include "Riostream.h"
#include "TEnv.h"
#include "TGListTree.h"
#include "TOrdCollection.h"
#include "TArrayF.h"
#include "TGHtml.h"
#include "TPRegexp.h"

#include "TEveManager.h"
#include "TEveViewer.h"
#include "TEveBrowser.h"
#include "TEveProjectionManager.h"
#include "TEveProjectionAxes.h"
#include "TEveScene.h"
#include "TEveGeoNode.h"
#include "TEveEventManager.h"
#include "TEveTrack.h"
#include "TEveSelection.h"

#include "TRootEmbeddedCanvas.h"
#include "TGSplitFrame.h"
#include "TGLOverlayButton.h"
#include "TGLEmbeddedViewer.h"
#include "TGDockableFrame.h"
#include "TGShapedFrame.h"
#include "TGButton.h"
#include "TGTab.h"

#include "TCanvas.h"
#include "TFormula.h"
#include "TF1.h"
#include "TH1F.h"

#ifdef WIN32
#include <TWin32SplashThread.h>
#endif

const char *filetypes[] = {
   "ROOT files",    "*.root",
   "All files",     "*",
   0,               0
};

const char *rcfiletypes[] = {
   "All files",     "*",
   0,               0
};

////////////////////////////////////////////////////////////////////////////////
class TGShapedToolTip : public TGShapedFrame {

private:
   TGShapedToolTip(const TGShapedToolTip&); // Not implemented
   TGShapedToolTip& operator=(const TGShapedToolTip&); // Not implemented

protected:
   Int_t                 fTextX, fTextY, fTextH;
   TString               fTextCol;

   TRootEmbeddedCanvas  *fEc;       // embedded canvas for histogram
   TH1                  *fHist;     // user histogram
   TString               fText;     // info (as tool tip) text

   virtual void          DoRedraw() {}

public:
   TGShapedToolTip(const char *picname, Int_t cx=0, Int_t cy=0, Int_t cw=0,
                   Int_t ch=0, Int_t tx=0, Int_t ty=0, Int_t th=0,
                   const char *col="#ffffff");
   virtual ~TGShapedToolTip();

   virtual void   CloseWindow();
   void           CreateCanvas(Int_t cx, Int_t cy, Int_t cw, Int_t ch);
   void           CreateCanvas(Int_t cw, Int_t ch, TGLayoutHints *hints);
   TH1           *GetHisto() const { return fHist; }
   const char    *GetText() const { return fText.Data(); }
   void           Refresh();
   void           SetHisto(TH1 *hist);
   void           SetText(const char *text);
   void           SetTextColor(const char *col);
   void           SetTextAttributes(Int_t tx, Int_t ty, Int_t th, const char *col=0);
   void           Show(Int_t x, Int_t y, const char *text = 0, TH1 *hist = 0);

   ClassDef(TGShapedToolTip, 0) // Shaped composite frame
};

////////////////////////////////////////////////////////////////////////////////
class HtmlObjTable : public TObject {
public:                     // make them public for shorter code

   TString   fName;
   Int_t     fNValues;      // number of values
   Int_t     fNFields;      // number of fields
   TArrayF  *fValues;
   TString  *fLabels;
   Bool_t    fExpand;

   TString   fHtml;         // HTML output code

   void Build();
   void BuildTitle();
   void BuildLabels();
   void BuildTable();

public:
   HtmlObjTable(const char *name, Int_t nfields, Int_t nvals, Bool_t exp=kTRUE);
   virtual ~HtmlObjTable();

   void     SetLabel(Int_t col, const char *label) { fLabels[col] = label; }
   void     SetValue(Int_t col, Int_t row, Float_t val) { fValues[col].SetAt(val, row); }
   TString  Html() const { return fHtml; }

   ClassDef(HtmlObjTable, 0);
};

////////////////////////////////////////////////////////////////////////////////
class HtmlSummary {
public:                           // make them public for shorter code
   Int_t           fNTables;
   TOrdCollection *fObjTables;    // ->array of object tables
   TString         fHtml;         // output HTML string
   TString         fTitle;        // page title
   TString         fHeader;       // HTML header
   TString         fFooter;       // HTML footer

   void     MakeHeader();
   void     MakeFooter();

public:
   HtmlSummary(const char *title);
   virtual ~HtmlSummary();

   HtmlObjTable  *AddTable(const char *name, Int_t nfields, Int_t nvals,
                           Bool_t exp=kTRUE, Option_t *opt="");
   HtmlObjTable  *GetTable(Int_t at) const { return (HtmlObjTable *)fObjTables->At(at); }
   void           Build();
   void           Clear(Option_t *option="");
   void           Reset(Option_t *option="");
   TString        Html() const { return fHtml; }

   ClassDef(HtmlSummary, 0);
};

////////////////////////////////////////////////////////////////////////////////
class SplitGLView : public TGMainFrame {

public:
   enum EMyCommands {
      kFileOpen, kFileExit, kFileLoadConfig, kFileSaveConfig,
      kHelpAbout, kGLPerspYOZ, kGLPerspXOZ, kGLPerspXOY, kGLXOY,
      kGLXOZ, kGLZOY, kGLOrthoRotate, kGLOrthoDolly, kSceneUpdate,
      kSceneUpdateAll, kSummaryUpdate
   };

private:
   TEvePad               *fPad;           // pad used as geometry container
   TGSplitFrame          *fSplitFrame;    // main (first) split frame
   TGLEmbeddedViewer     *fViewer0;       // main GL viewer
   TGLEmbeddedViewer     *fViewer1;       // first GL viewer
   TGLEmbeddedViewer     *fViewer2;       // second GL viewer
   TGLEmbeddedViewer     *fActViewer;     // actual (active) GL viewer
   static HtmlSummary    *fgHtmlSummary;  // summary HTML table
   static TGHtml         *fgHtml;
   TGMenuBar             *fMenuBar;       // main menu bar
   TGPopupMenu           *fMenuFile;      // 'File' popup menu
   TGPopupMenu           *fMenuHelp;      // 'Help' popup menu
   TGPopupMenu           *fMenuCamera;    // 'Camera' popup menu
   TGPopupMenu           *fMenuScene;     // 'Scene' popup menu
   TGStatusBar           *fStatusBar;     // status bar
   TGShapedToolTip       *fShapedToolTip; // shaped tooltip
   Bool_t                 fIsEmbedded;

   TEveViewer            *fViewer[3];
   TEveProjectionManager *fRPhiMgr;
   TEveProjectionManager *fRhoZMgr;

public:
   SplitGLView(const TGWindow *p=0, UInt_t w=800, UInt_t h=600, Bool_t embed=kFALSE);
   virtual ~SplitGLView();

   void           ItemClicked(TGListTreeItem *item, Int_t btn, Int_t x, Int_t y);
   void           HandleMenu(Int_t id);
   void           OnClicked(TObject *obj);
   void           OnMouseIdle(TGLPhysicalShape *shape, UInt_t posx, UInt_t posy);
   void           OnMouseOver(TGLPhysicalShape *shape);
   void           OnViewerActivated();
   void           OpenFile(const char *fname);
   void           SwapToMainView(TGLViewerBase *viewer);
   void           ToggleOrthoRotate();
   void           ToggleOrthoDolly();
   void           UnDock(TGLViewerBase *viewer);
   void           LoadConfig(const char *fname);
   void           SaveConfig(const char *fname);
   static void    UpdateSummary();

   TEveProjectionManager *GetRPhiMgr() const { return fRPhiMgr; }
   TEveProjectionManager *GetRhoZMgr() const { return fRhoZMgr; }

   ClassDef(SplitGLView, 0)
};

TEveProjectionManager *gRPhiMgr = 0;
TEveProjectionManager *gRhoZMgr = 0;

ClassImp(TGShapedToolTip)
ClassImp(HtmlObjTable)
ClassImp(HtmlSummary)
ClassImp(SplitGLView)

HtmlSummary *SplitGLView::fgHtmlSummary = 0;
TGHtml *SplitGLView::fgHtml = 0;

//______________________________________________________________________________
TGShapedToolTip::TGShapedToolTip(const char *pname, Int_t cx, Int_t cy, Int_t cw,
                             Int_t ch, Int_t tx, Int_t ty, Int_t th,
                             const char *col) :
   TGShapedFrame(pname, gClient->GetDefaultRoot(), 400, 300, kTempFrame |
                 kHorizontalFrame), fEc(0), fHist(0)
{
   // Shaped window constructor

   fTextX = tx; fTextY = ty; fTextH = th;
   if (col)
      fTextCol = col;
   else
      fTextCol = "0x000000";

   // create the embedded canvas
   if ((cx > 0) && (cy > 0) && (cw > 0) && (ch > 0)) {
      Int_t lhRight  = fWidth-cx-cw;
      Int_t lhBottom = fHeight-cy-ch;
      fEc = new TRootEmbeddedCanvas("ec", this, cw, ch, 0);
      AddFrame(fEc, new TGLayoutHints(kLHintsTop | kLHintsLeft, cx,
                                      lhRight, cy, lhBottom));
   }
   MapSubwindows();
   Resize();
   Resize(fBgnd->GetWidth(), fBgnd->GetHeight());
}

//______________________________________________________________________________
TGShapedToolTip::~TGShapedToolTip()
{
   // Destructor.

   if (fHist)
      delete fHist;
   if (fEc)
      delete fEc;
}

//______________________________________________________________________________
void TGShapedToolTip::CloseWindow()
{
   // Close shaped window.

   DeleteWindow();
}

//______________________________________________________________________________
void TGShapedToolTip::Refresh()
{
   // Redraw the window with current attributes.

   const char *str = fText.Data();
   char *string = strdup(str);
   Int_t nlines = 0, size = fTextH;
   TString ar = "arial.ttf";
   char *s = strtok((char *)string, "\n");
   TImage *img = (TImage*)fImage->Clone("img");
   img->DrawText(fTextX, fTextY+(nlines*size), s, size, fTextCol, ar);
   while ((s = strtok(0, "\n"))) {
      nlines++;
      img->DrawText(fTextX, fTextY+(nlines*size), s, size, fTextCol, ar);
   }
   img->PaintImage(fId, 0, 0, 0, 0, 0, 0, "opaque");
   free(string);
   delete img;
   gVirtualX->Update();
}

//______________________________________________________________________________
void TGShapedToolTip::CreateCanvas(Int_t cx, Int_t cy, Int_t cw, Int_t ch)
{

   // create the embedded canvas
   Int_t lhRight  = fWidth-cx-cw;
   Int_t lhBottom = fHeight-cy-ch;
   fEc = new TRootEmbeddedCanvas("ec", this, cw, ch, 0);
   AddFrame(fEc, new TGLayoutHints(kLHintsTop | kLHintsLeft, cx,
                                   lhRight, cy, lhBottom));
   MapSubwindows();
   Resize();
   Resize(fBgnd->GetWidth(), fBgnd->GetHeight());
   if (IsMapped()) {
      Refresh();
   }
}

//______________________________________________________________________________
void TGShapedToolTip::CreateCanvas(Int_t cw, Int_t ch, TGLayoutHints *hints)
{
   // Create the embedded canvas.

   fEc = new TRootEmbeddedCanvas("ec", this, cw, ch, 0);
   AddFrame(fEc, hints);
   MapSubwindows();
   Resize();
   Resize(fBgnd->GetWidth(), fBgnd->GetHeight());
   if (IsMapped()) {
      Refresh();
   }
}

//______________________________________________________________________________
void TGShapedToolTip::SetHisto(TH1 *hist)
{
   // Set which histogram has to be displayed in the embedded canvas.

   if (hist) {
      if (fHist) {
         delete fHist;
         if (fEc)
            fEc->GetCanvas()->Clear();
      }
      fHist = (TH1 *)hist->Clone();
      if (fEc) {
         fEc->GetCanvas()->SetBorderMode(0);
         fEc->GetCanvas()->SetFillColor(10);
         fEc->GetCanvas()->cd();
         fHist->Draw();
         fEc->GetCanvas()->Update();
      }
   }
}

//______________________________________________________________________________
void TGShapedToolTip::SetText(const char *text)
{
   // Set which text has to be displayed.

   if (text) {
      fText = text;
   }
   if (IsMapped())
      Refresh();
}

//______________________________________________________________________________
void TGShapedToolTip::SetTextColor(const char *col)
{
   // Set text color.

   fTextCol = col;
   if (IsMapped())
      Refresh();
}

//______________________________________________________________________________
void TGShapedToolTip::SetTextAttributes(Int_t tx, Int_t ty, Int_t th,
                                        const char *col)
{
   // Set text attributes (position, size and color).

   fTextX = tx; fTextY = ty; fTextH = th;
   if (col)
      fTextCol = col;
   if (IsMapped())
      Refresh();
}

//______________________________________________________________________________
void TGShapedToolTip::Show(Int_t x, Int_t y, const char *text, TH1 *hist)
{
   // Show (popup) the shaped window at location x,y and possibly
   // set the text and histogram to be displayed.

   Move(x, y);
   MapWindow();

   if (text)
      SetText(text);
   if (hist)
      SetHisto(hist);
   // end of demo code -------------------------------------------
   if (fHist) {
      fEc->GetCanvas()->SetBorderMode(0);
      fEc->GetCanvas()->SetFillColor(10);
      fEc->GetCanvas()->cd();
      fHist->Draw();
      fEc->GetCanvas()->Update();
   }
   Refresh();
}

//______________________________________________________________________________
HtmlObjTable::HtmlObjTable(const char *name, Int_t nfields, Int_t nvals, Bool_t exp) :
   fName(name), fNValues(nvals), fNFields(nfields), fExpand(exp)
{
   // Constructor.

   fValues = new TArrayF[fNFields];
   for (int i=0;i<fNFields;i++)
      fValues[i].Set(nvals);
   fLabels = new TString[fNFields];
}

//______________________________________________________________________________
HtmlObjTable::~HtmlObjTable()
{
   // Destructor.

   delete [] fValues;
   delete [] fLabels;
}

//______________________________________________________________________________
void HtmlObjTable::Build()
{
   // Build HTML code.

   fHtml = "<table width=100% border=1 cellspacing=0 cellpadding=0 bgcolor=f0f0f0> ",

   BuildTitle();
   if (fExpand && (fNFields > 0) && (fNValues > 0)) {
      BuildLabels();
      BuildTable();
   }

   fHtml += "</table>";
}

//______________________________________________________________________________
void HtmlObjTable::BuildTitle()
{
   // Build table title.

   fHtml += "<tr><td colspan=";
   fHtml += Form("%d>", fNFields+1);
   fHtml += "<table width=100% border=0 cellspacing=2 cellpadding=0 bgcolor=6e6ea0>";
   fHtml += "<tr><td align=left>";
   fHtml += "<font face=Verdana size=3 color=ffffff><b><i>";
   fHtml += fName;
   fHtml += "</i></b></font></td>";
   fHtml += "<td>";
   fHtml += "<td align=right> ";
   fHtml += "<font face=Verdana size=3 color=ffffff><b><i>";
   fHtml += Form("Size = %d", fNValues);
   fHtml += "</i></b></font></td></tr>";
   fHtml += "</table>";
   fHtml += "</td></tr>";
}

//______________________________________________________________________________
void HtmlObjTable::BuildLabels()
{
   // Build table labels.

   Int_t i;
   fHtml += "<tr bgcolor=c0c0ff>";
   fHtml += "<th> </th>"; // for the check boxes
   for (i=0;i<fNFields;i++) {
      fHtml += "<th> ";
      fHtml += fLabels[i];
      fHtml += " </th>"; // for the check boxes
   }
   fHtml += "</tr>";
}

//______________________________________________________________________________
void HtmlObjTable::BuildTable()
{
   // Build part of table with values.

   for (int i = 0; i < fNValues; i++) {
      if (i%2)
         fHtml += "<tr bgcolor=e0e0ff>";
      else
         fHtml += "<tr bgcolor=ffffff>";

      TString name = fName;
      name.ReplaceAll(" ", "_");
      // checkboxes
      fHtml += "<td bgcolor=d0d0ff align=\"center\">";
      fHtml += "<input type=\"checkbox\" name=\"";
      fHtml += name;
      fHtml += Form("[%d]\">",i);
      fHtml += "</td>";

      for (int j = 0; j < fNFields; j++) {
         fHtml += "<td width=";
         fHtml += Form("%d%%", 100/fNFields);
         fHtml += " align=\"center\"";
         fHtml += ">";
         fHtml += Form("%1.4f", fValues[j][i]);
         fHtml += "</td>";
      }
      fHtml += "</tr> ";
   }
}

//______________________________________________________________________________
HtmlSummary::HtmlSummary(const char *title) : fNTables(0), fTitle(title)
{
   // Constructor.

   fObjTables = new TOrdCollection();
}

//______________________________________________________________________________
HtmlSummary::~HtmlSummary()
{
   // Destructor.

   Reset();
}

//______________________________________________________________________________
HtmlObjTable *HtmlSummary::AddTable(const char *name, Int_t nfields, Int_t nvals,
                                    Bool_t exp, Option_t *option)
{
   // Add a new table in our list of tables.

   TString opt = option;
   opt.ToLower();
   HtmlObjTable *table = new HtmlObjTable(name, nfields, nvals, exp);
   fNTables++;
   if (opt.Contains("first"))
      fObjTables->AddFirst(table);
   else
      fObjTables->Add(table);
   return table;
}

//______________________________________________________________________________
void HtmlSummary::Clear(Option_t *option)
{
   // Clear the table list.

   if (option && option[0] == 'D')
      fObjTables->Delete(option);
   else
      fObjTables->Clear(option);
   fNTables = 0;
}

//______________________________________________________________________________
void HtmlSummary::Reset(Option_t *)
{
   // Reset (delete) the table list;

   delete fObjTables; fObjTables = 0;
   fNTables = 0;
}

//______________________________________________________________________________
void HtmlSummary::Build()
{
   // Build the summary.

   MakeHeader();
   for (int i=0;i<fNTables;i++) {
      GetTable(i)->Build();
      fHtml += GetTable(i)->Html();
   }
   MakeFooter();
}

//______________________________________________________________________________
void HtmlSummary::MakeHeader()
{
   // Make HTML header.

   fHeader  = "<html><head><title>";
   fHeader += fTitle;
   fHeader += "</title></head><body>";
   fHeader += "<center><h2><font color=#2222ee><i>";
   fHeader += fTitle;
   fHeader += "</i></font></h2></center>";
   fHtml    = fHeader;
}

//______________________________________________________________________________
void HtmlSummary::MakeFooter()
{
   // Make HTML footer.

   fFooter  = "<br><p><br><center><strong><font size=2 color=#2222ee>";
   fFooter += "Example of using Html widget to display tabular data";
   fFooter += "<br>";
   fFooter += "(c) 2007-2010 Bertrand Bellenot";
   fFooter += "</font></strong></center></body></html>";
   fHtml   += fFooter;
}

//______________________________________________________________________________
SplitGLView::SplitGLView(const TGWindow *p, UInt_t w, UInt_t h, Bool_t embed) :
   TGMainFrame(p, w, h), fActViewer(0), fShapedToolTip(0), fIsEmbedded(embed)
{
   // Main frame constructor.

   TGSplitFrame *frm;
   TEveScene *s = 0;
   TGHorizontalFrame *hfrm;
   TGDockableFrame *dfrm;
   TGPictureButton *button;

   // create the "file" popup menu
   fMenuFile = new TGPopupMenu(gClient->GetRoot());
   fMenuFile->AddEntry("&Open...", kFileOpen);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry( "&Update Summary", kSummaryUpdate);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("&Load Config...", kFileLoadConfig);
   fMenuFile->AddEntry("&Save Config...", kFileSaveConfig);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("E&xit", kFileExit);

   // create the "camera" popup menu
   fMenuCamera = new TGPopupMenu(gClient->GetRoot());
   fMenuCamera->AddEntry("Perspective (Floor XOZ)", kGLPerspXOZ);
   fMenuCamera->AddEntry("Perspective (Floor YOZ)", kGLPerspYOZ);
   fMenuCamera->AddEntry("Perspective (Floor XOY)", kGLPerspXOY);
   fMenuCamera->AddEntry("Orthographic (XOY)", kGLXOY);
   fMenuCamera->AddEntry("Orthographic (XOZ)", kGLXOZ);
   fMenuCamera->AddEntry("Orthographic (ZOY)", kGLZOY);
   fMenuCamera->AddSeparator();
   fMenuCamera->AddEntry("Ortho allow rotate", kGLOrthoRotate);
   fMenuCamera->AddEntry("Ortho allow dolly",  kGLOrthoDolly);

   fMenuScene = new TGPopupMenu(gClient->GetRoot());
   fMenuScene->AddEntry("&Update Current", kSceneUpdate);
   fMenuScene->AddEntry("Update &All", kSceneUpdateAll);

   // create the "help" popup menu
   fMenuHelp = new TGPopupMenu(gClient->GetRoot());
   fMenuHelp->AddEntry("&About", kHelpAbout);

   // create the main menu bar
   fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
   fMenuBar->AddPopup("&File", fMenuFile, new TGLayoutHints(kLHintsTop |
                      kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Camera", fMenuCamera, new TGLayoutHints(kLHintsTop |
                      kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Scene", fMenuScene, new TGLayoutHints(kLHintsTop |
                      kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Help", fMenuHelp, new TGLayoutHints(kLHintsTop |
                      kLHintsRight));

   AddFrame(fMenuBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   // connect menu signals to our menu handler slot
   fMenuFile->Connect("Activated(Int_t)", "SplitGLView", this,
                      "HandleMenu(Int_t)");
   fMenuCamera->Connect("Activated(Int_t)", "SplitGLView", this,
                        "HandleMenu(Int_t)");
   fMenuScene->Connect("Activated(Int_t)", "SplitGLView", this,
                       "HandleMenu(Int_t)");
   fMenuHelp->Connect("Activated(Int_t)", "SplitGLView", this,
                      "HandleMenu(Int_t)");

   if (fIsEmbedded && gEve) {
      // use status bar from the browser
      fStatusBar = gEve->GetBrowser()->GetStatusBar();
   }
   else {
      // create the status bar
      Int_t parts[] = {45, 15, 10, 30};
      fStatusBar = new TGStatusBar(this, 50, 10);
      fStatusBar->SetParts(parts, 4);
      AddFrame(fStatusBar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX,
               0, 0, 10, 0));
   }

   // create eve pad (our geometry container)
   fPad = new TEvePad();
   fPad->SetFillColor(kBlack);

   // create the split frames
   fSplitFrame = new TGSplitFrame(this, 800, 600);
   AddFrame(fSplitFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   // split it once
   fSplitFrame->HSplit(434);
   // then split each part again (this will make four parts)
   fSplitFrame->GetSecond()->VSplit(266);
   fSplitFrame->GetSecond()->GetSecond()->VSplit(266);

   TGLOverlayButton *but1, *but2, *but3, *but4, *but5, *but6;
   // get top (main) split frame
   frm = fSplitFrame->GetFirst();
   frm->SetName("Main_View");

   // create (embed) a GL viewer inside
   fViewer0 = new TGLEmbeddedViewer(frm, fPad);
   but1 = new TGLOverlayButton(fViewer0, "Swap", 10.0, -10.0, 55.0, 16.0);
   but1->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "SwapToMainView(TGLViewerBase*)");
   but2 = new TGLOverlayButton(fViewer0, "Undock", 70.0, -10.0, 55.0, 16.0);
   but2->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "UnDock(TGLViewerBase*)");
   frm->AddFrame(fViewer0->GetFrame(), new TGLayoutHints(kLHintsExpandX |
                 kLHintsExpandY));
   // set the camera to perspective (XOZ) for this viewer
   fViewer0->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
   // connect signal we are interested to
   fViewer0->Connect("MouseOver(TGLPhysicalShape*)", "SplitGLView", this,
                      "OnMouseOver(TGLPhysicalShape*)");
   fViewer0->Connect("Activated()", "SplitGLView", this,
                      "OnViewerActivated()");
   fViewer0->Connect("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)",
                      "SplitGLView", this,
                      "OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");
   fViewer0->Connect("Clicked(TObject*)", "SplitGLView", this,
                      "OnClicked(TObject*)");
   fViewer[0] = new TEveViewer("SplitGLViewer[0]");
   fViewer[0]->SetGLViewer(fViewer0, fViewer0->GetFrame());
   fViewer[0]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
      fViewer[0]->AddScene(gEve->GetGlobalScene());
      fViewer[0]->AddScene(gEve->GetEventScene());
      gEve->GetViewers()->AddElement(fViewer[0]);
      s = gEve->SpawnNewScene("Rho-Z Projection");
      // projections
      fRhoZMgr = new TEveProjectionManager(TEveProjection::kPT_RhoZ);
      s->AddElement(fRhoZMgr);
      gEve->AddToListTree(fRhoZMgr, kTRUE);
      TEveProjectionAxes* a = new TEveProjectionAxes(fRhoZMgr);
      s->AddElement(a);
   }

   // get bottom left split frame
   frm = fSplitFrame->GetSecond()->GetFirst();
   frm->SetName("Bottom_Left");

   // create (embed) a GL viewer inside
   fViewer1 = new TGLEmbeddedViewer(frm, fPad);
   but3 = new TGLOverlayButton(fViewer1, "Swap", 10.0, -10.0, 55.0, 16.0);
   but3->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "SwapToMainView(TGLViewerBase*)");
   but4 = new TGLOverlayButton(fViewer1, "Undock", 70.0, -10.0, 55.0, 16.0);
   but4->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "UnDock(TGLViewerBase*)");
   frm->AddFrame(fViewer1->GetFrame(), new TGLayoutHints(kLHintsExpandX |
                  kLHintsExpandY));

   // set the camera to orthographic (XOY) for this viewer
   fViewer1->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   // connect signal we are interested to
   fViewer1->Connect("MouseOver(TGLPhysicalShape*)", "SplitGLView", this,
                      "OnMouseOver(TGLPhysicalShape*)");
   fViewer1->Connect("Activated()", "SplitGLView", this,
                      "OnViewerActivated()");
   fViewer1->Connect("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)",
                      "SplitGLView", this,
                      "OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");
   fViewer1->Connect("Clicked(TObject*)", "SplitGLView", this,
                      "OnClicked(TObject*)");
   fViewer[1] = new TEveViewer("SplitGLViewer[1]");
   fViewer[1]->SetGLViewer(fViewer1, fViewer1->GetFrame());
   fViewer[1]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
      fRhoZMgr->ImportElements((TEveElement *)gEve->GetGlobalScene());
      fRhoZMgr->ImportElements((TEveElement *)gEve->GetEventScene());
      fViewer[1]->AddScene(s);
      gEve->GetViewers()->AddElement(fViewer[1]);
      gRhoZMgr = fRhoZMgr;

      s = gEve->SpawnNewScene("R-Phi Projection");
      // projections
      fRPhiMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
      s->AddElement(fRPhiMgr);
      gEve->AddToListTree(fRPhiMgr, kTRUE);
      TEveProjectionAxes* a = new TEveProjectionAxes(fRPhiMgr);
      s->AddElement(a);
   }

   // get bottom center split frame
   frm = fSplitFrame->GetSecond()->GetSecond()->GetFirst();
   frm->SetName("Bottom_Center");

   // create (embed) a GL viewer inside
   fViewer2 = new TGLEmbeddedViewer(frm, fPad);
   but5 = new TGLOverlayButton(fViewer2, "Swap", 10.0, -10.0, 55.0, 16.0);
   but5->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "SwapToMainView(TGLViewerBase*)");
   but6 = new TGLOverlayButton(fViewer2, "Undock", 70.0, -10.0, 55.0, 16.0);
   but6->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "UnDock(TGLViewerBase*)");
   frm->AddFrame(fViewer2->GetFrame(), new TGLayoutHints(kLHintsExpandX |
                  kLHintsExpandY));

   // set the camera to orthographic (XOY) for this viewer
   fViewer2->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   // connect signal we are interested to
   fViewer2->Connect("MouseOver(TGLPhysicalShape*)", "SplitGLView", this,
                      "OnMouseOver(TGLPhysicalShape*)");
   fViewer2->Connect("Activated()", "SplitGLView", this,
                      "OnViewerActivated()");
   fViewer2->Connect("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)",
                      "SplitGLView", this,
                      "OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");
   fViewer2->Connect("Clicked(TObject*)", "SplitGLView", this,
                      "OnClicked(TObject*)");
   fViewer[2] = new TEveViewer("SplitGLViewer[2]");
   fViewer[2]->SetGLViewer(fViewer2, fViewer2->GetFrame());
   fViewer[2]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
      fRPhiMgr->ImportElements((TEveElement *)gEve->GetGlobalScene());
      fRPhiMgr->ImportElements((TEveElement *)gEve->GetEventScene());
      fViewer[2]->AddScene(s);
      gEve->GetViewers()->AddElement(fViewer[2]);
      gRPhiMgr = fRPhiMgr;
   }

   // get bottom right split frame
   frm = fSplitFrame->GetSecond()->GetSecond()->GetSecond();
   frm->SetName("Bottom_Right");

   dfrm = new TGDockableFrame(frm);
   dfrm->SetFixedSize(kFALSE);
   dfrm->EnableHide(kFALSE);
   hfrm = new TGHorizontalFrame(dfrm);
   button= new TGPictureButton(hfrm, gClient->GetPicture("swap.png"));
   button->SetToolTipText("Swap to big view");
   hfrm->AddFrame(button);
   button->Connect("Clicked()","SplitGLView",this,"SwapToMainView(TGLViewerBase*=0)");
   fgHtmlSummary = new HtmlSummary("Alice Event Display Summary Table");
   fgHtml = new TGHtml(hfrm, 100, 100, -1);
   hfrm->AddFrame(fgHtml, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   dfrm->AddFrame(hfrm, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   frm->AddFrame(dfrm, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   if (fIsEmbedded && gEve) {
      gEve->GetListTree()->Connect("Clicked(TGListTreeItem*, Int_t, Int_t, Int_t)",
         "SplitGLView", this, "ItemClicked(TGListTreeItem*, Int_t, Int_t, Int_t)");
   }

   fShapedToolTip = new TGShapedToolTip("Default.png", 120, 22, 160, 110,
                                        23, 115, 12, "#ffff80");
   Resize(GetDefaultSize());
   MapSubwindows();
   MapWindow();
   LoadConfig(".everc");
}

//______________________________________________________________________________
SplitGLView::~SplitGLView()
{
   // Clean up main frame...
   //Cleanup();

   fMenuFile->Disconnect("Activated(Int_t)", this, "HandleMenu(Int_t)");
   fMenuCamera->Disconnect("Activated(Int_t)", this, "HandleMenu(Int_t)");
   fMenuScene->Disconnect("Activated(Int_t)", this, "HandleMenu(Int_t)");
   fMenuHelp->Disconnect("Activated(Int_t)", this, "HandleMenu(Int_t)");
   fViewer0->Disconnect("MouseOver(TGLPhysicalShape*)", this,
                         "OnMouseOver(TGLPhysicalShape*)");
   fViewer0->Disconnect("Activated()", this, "OnViewerActivated()");
   fViewer0->Disconnect("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)",
                         this, "OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");
   fViewer1->Disconnect("MouseOver(TGLPhysicalShape*)", this,
                         "OnMouseOver(TGLPhysicalShape*)");
   fViewer1->Disconnect("Activated()", this, "OnViewerActivated()");
   fViewer1->Disconnect("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)",
                         this, "OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");
   fViewer2->Disconnect("MouseOver(TGLPhysicalShape*)", this,
                         "OnMouseOver(TGLPhysicalShape*)");
   fViewer2->Disconnect("Activated()", this, "OnViewerActivated()");
   fViewer2->Disconnect("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)",
                         this, "OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");
   if (!fIsEmbedded) {
      delete fViewer[0];
      delete fViewer[1];
      delete fViewer[2];
   }
   delete fShapedToolTip;
   delete fMenuFile;
   delete fMenuScene;
   delete fMenuCamera;
   delete fMenuHelp;
   if (!fIsEmbedded)
      delete fMenuBar;
   delete fViewer0;
   delete fViewer1;
   delete fViewer2;
   delete fSplitFrame;
   delete fPad;
   if (!fIsEmbedded) {
      delete fStatusBar;
      gApplication->Terminate(0);
   }
}

//______________________________________________________________________________
void SplitGLView::HandleMenu(Int_t id)
{
   // Handle menu items.

   static TString rcdir(".");
   static TString rcfile(".everc");

   switch (id) {

      case kFileOpen:
         {
            static TString dir(".");
            TGFileInfo fi;
            fi.fFileTypes = filetypes;
            fi.SetIniDir(dir);
            new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
            if (fi.fFilename)
               OpenFile(fi.fFilename);
            dir = fi.fIniDir;
         }
         break;

      case kFileLoadConfig:
         {
            TGFileInfo fi;
            fi.fFileTypes = rcfiletypes;
            fi.SetIniDir(rcdir);
            fi.SetFilename(rcfile);
            new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
            if (fi.fFilename) {
               rcfile = fi.fFilename;
               LoadConfig(fi.fFilename);
            }
            rcdir = fi.fIniDir;
         }
         break;

      case kFileSaveConfig:
         {
            TGFileInfo fi;
            fi.fFileTypes = rcfiletypes;
            fi.SetIniDir(rcdir);
            fi.SetFilename(rcfile);
            new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fi);
            if (fi.fFilename) {
               rcfile = fi.fFilename;
               SaveConfig(fi.fFilename);
            }
            rcdir = fi.fIniDir;
         }
         break;

      case kFileExit:
         CloseWindow();
         break;

      case kGLPerspYOZ:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraPerspYOZ);
         break;
      case kGLPerspXOZ:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
         break;
      case kGLPerspXOY:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraPerspXOY);
         break;
      case kGLXOY:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
         break;
      case kGLXOZ:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
         break;
      case kGLZOY:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
         break;
      case kGLOrthoRotate:
         ToggleOrthoRotate();
         break;
      case kGLOrthoDolly:
         ToggleOrthoDolly();
         break;

      case kSceneUpdate:
         if (fActViewer)
            fActViewer->UpdateScene();
         UpdateSummary();
         break;

      case kSceneUpdateAll:
         fViewer0->UpdateScene();
         fViewer1->UpdateScene();
         fViewer2->UpdateScene();
         UpdateSummary();
         break;

      case kSummaryUpdate:
         UpdateSummary();
         break;

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

      default:
         break;
   }
}

//______________________________________________________________________________
void SplitGLView::OnClicked(TObject *obj)
{
   // Handle click events in GL viewer

   if (obj)
      fStatusBar->SetText(Form("User clicked on: \"%s\"", obj->GetName()), 1);
   else
      fStatusBar->SetText("", 1);
}

//______________________________________________________________________________
void SplitGLView::OnMouseIdle(TGLPhysicalShape *shape, UInt_t posx, UInt_t posy)
{
   // Slot used to handle "OnMouseIdle" signal coming from any GL viewer.
   // We receive a pointer on the physical shape in which the mouse cursor is
   // and the actual cursor position (x,y)

   Window_t wtarget;
   Int_t    x = 0, y = 0;

   static TH1F *h1f = 0;
   TFormula *form1 = new TFormula("form1","abs(sin(x)/x)");
   TF1 *sqroot = new TF1("sqroot","x*gaus(0) + [3]*form1",0,10);
   sqroot->SetParameters(10,4,1,20);
   if (h1f == 0)
      h1f = new TH1F("h1f","",50,0,10);
   h1f->Reset();
   h1f->SetFillColor(45);
   h1f->SetStats(0);
   h1f->FillRandom("sqroot",200);

   if (fShapedToolTip) {
      fShapedToolTip->UnmapWindow();
   }
   if (shape && shape->GetLogical() && shape->GetLogical()->GetExternal()) {
      // get the actual viewer who actually emitted the signal
      TGLEmbeddedViewer *actViewer = dynamic_cast<TGLEmbeddedViewer*>((TQObject*)gTQSender);
      // then translate coordinates from the root (screen) coordinates
      // to the actual frame (viewer) ones
      gVirtualX->TranslateCoordinates(actViewer->GetFrame()->GetId(),
               gClient->GetDefaultRoot()->GetId(), posx, posy, x, y,
               wtarget);
      // Then display our tooltip at this x,y location
      if (fShapedToolTip) {
         fShapedToolTip->Show(x+5, y+5, Form("%s\n     \n%s",
                              shape->GetLogical()->GetExternal()->IsA()->GetName(),
                              shape->GetLogical()->GetExternal()->GetName()), h1f);
      }
   }
}

//______________________________________________________________________________
void SplitGLView::OnMouseOver(TGLPhysicalShape *shape)
{
   // Slot used to handle "OnMouseOver" signal coming from any GL viewer.
   // We receive a pointer on the physical shape in which the mouse cursor is.

   // display information on the physical shape in the status bar
   if (shape && shape->GetLogical() && shape->GetLogical()->GetExternal())
      fStatusBar->SetText(Form("Mouse Over: \"%s\"",
         shape->GetLogical()->GetExternal()->GetName()), 0);
   else
      fStatusBar->SetText("", 0);
}

//______________________________________________________________________________
void SplitGLView::OnViewerActivated()
{
   // Slot used to handle "Activated" signal coming from any GL viewer.
   // Used to know which GL viewer is active.

   static Pixel_t green = 0;
   // set the actual GL viewer frame to default color
   if (fActViewer && fActViewer->GetFrame())
      fActViewer->GetFrame()->ChangeBackground(GetDefaultFrameBackground());

   // change the actual GL viewer to the one who emitted the signal
   // fActViewer = (TGLEmbeddedViewer *)gTQSender;
   fActViewer = dynamic_cast<TGLEmbeddedViewer*>((TQObject*)gTQSender);

   if (fActViewer == 0) {
      printf ("dyncast failed ...\n");
      return;
   }

   // get the highlight color (only once)
   if (green == 0) {
      gClient->GetColorByName("green", green);
   }
   // set the new actual GL viewer frame to highlight color
   if (fActViewer->GetFrame())
      fActViewer->GetFrame()->ChangeBackground(green);

   // update menu entries to match actual viewer's options
   if (fActViewer->GetOrthoXOYCamera()->GetDollyToZoom() &&
       fActViewer->GetOrthoXOZCamera()->GetDollyToZoom() &&
       fActViewer->GetOrthoZOYCamera()->GetDollyToZoom())
      fMenuCamera->UnCheckEntry(kGLOrthoDolly);
   else
      fMenuCamera->CheckEntry(kGLOrthoDolly);

   if (fActViewer->GetOrthoXOYCamera()->GetEnableRotate() &&
       fActViewer->GetOrthoXOYCamera()->GetEnableRotate() &&
       fActViewer->GetOrthoXOYCamera()->GetEnableRotate())
      fMenuCamera->CheckEntry(kGLOrthoRotate);
   else
      fMenuCamera->UnCheckEntry(kGLOrthoRotate);
}

//______________________________________________________________________________
void SplitGLView::OpenFile(const char *fname)
{
   // Open a Root file to display a geometry in the GL viewers.

   TString filename = fname;
   // check if the file type is correct
   if (!filename.EndsWith(".root")) {
      new TGMsgBox(gClient->GetRoot(), this, "OpenFile",
                   Form("The file \"%s\" is not a root file!", fname),
                   kMBIconExclamation, kMBOk);
      return;
   }
   // check if the root file contains a geometry
   if (TGeoManager::Import(fname) == 0) {
      new TGMsgBox(gClient->GetRoot(), this, "OpenFile",
                   Form("The file \"%s\" does't contain a geometry", fname),
                   kMBIconExclamation, kMBOk);
      return;
   }
   gGeoManager->DefaultColors();
   // delete previous primitives (if any)
   fPad->GetListOfPrimitives()->Delete();
   // and add the geometry to eve pad (container)
   fPad->GetListOfPrimitives()->Add(gGeoManager->GetTopVolume());
   // paint the geometry in each GL viewer
   fViewer0->PadPaint(fPad);
   fViewer1->PadPaint(fPad);
   fViewer2->PadPaint(fPad);
}

//______________________________________________________________________________
void SplitGLView::ToggleOrthoRotate()
{
   // Toggle state of the 'Ortho allow rotate' menu entry.

   if (fMenuCamera->IsEntryChecked(kGLOrthoRotate))
      fMenuCamera->UnCheckEntry(kGLOrthoRotate);
   else
      fMenuCamera->CheckEntry(kGLOrthoRotate);
   Bool_t state = fMenuCamera->IsEntryChecked(kGLOrthoRotate);
   if (fActViewer) {
      fActViewer->GetOrthoXOYCamera()->SetEnableRotate(state);
      fActViewer->GetOrthoXOYCamera()->SetEnableRotate(state);
      fActViewer->GetOrthoXOYCamera()->SetEnableRotate(state);
   }
}

//______________________________________________________________________________
void SplitGLView::ToggleOrthoDolly()
{
   // Toggle state of the 'Ortho allow dolly' menu entry.

   if (fMenuCamera->IsEntryChecked(kGLOrthoDolly))
      fMenuCamera->UnCheckEntry(kGLOrthoDolly);
   else
      fMenuCamera->CheckEntry(kGLOrthoDolly);
   Bool_t state = ! fMenuCamera->IsEntryChecked(kGLOrthoDolly);
   if (fActViewer) {
      fActViewer->GetOrthoXOYCamera()->SetDollyToZoom(state);
      fActViewer->GetOrthoXOZCamera()->SetDollyToZoom(state);
      fActViewer->GetOrthoZOYCamera()->SetDollyToZoom(state);
   }
}

//______________________________________________________________________________
void SplitGLView::ItemClicked(TGListTreeItem *item, Int_t, Int_t, Int_t)
{
   // Item has been clicked, based on mouse button do:

   static const TEveException eh("SplitGLView::ItemClicked ");
   TEveElement* re = (TEveElement*)item->GetUserData();
   if(re == 0) return;
   TObject* obj = re->GetObject(eh);
   if (obj->InheritsFrom("TEveViewer")) {
      TGLViewer *v = ((TEveViewer *)obj)->GetGLViewer();
      //v->Activated();
      if (v->InheritsFrom("TGLEmbeddedViewer")) {
         TGLEmbeddedViewer *ev = (TGLEmbeddedViewer *)v;
         gVirtualX->SetInputFocus(ev->GetGLWidget()->GetId());
      }
   }
}

//______________________________________________________________________________
void SplitGLView::LoadConfig(const char *fname)
{

   Int_t height, width;
   TEnv *env = new TEnv(fname);

   Int_t mainheight = env->GetValue("MainView.Height", 434);
   Int_t blwidth    = env->GetValue("Bottom.Left.Width", 266);
   Int_t bcwidth    = env->GetValue("Bottom.Center.Width", 266);
   Int_t brwidth    = env->GetValue("Bottom.Right.Width", 266);
   Int_t top_height = env->GetValue("Right.Tab.Height", 0);
   Int_t bottom_height = env->GetValue("Bottom.Tab.Height", 0);

   if (fIsEmbedded && gEve) {
      Int_t sel = env->GetValue("Eve.Selection", gEve->GetSelection()->GetPickToSelect());
      Int_t hi = env->GetValue("Eve.Highlight", gEve->GetHighlight()->GetPickToSelect());
      gEve->GetBrowser()->EveMenu(9+sel);
      gEve->GetBrowser()->EveMenu(13+hi);

      width  = env->GetValue("Eve.Width", (Int_t)gEve->GetBrowser()->GetWidth());
      height = env->GetValue("Eve.Height", (Int_t)gEve->GetBrowser()->GetHeight());
      gEve->GetBrowser()->Resize(width, height);
   }

   // top (main) split frame
   width = fSplitFrame->GetFirst()->GetWidth();
   fSplitFrame->GetFirst()->Resize(width, mainheight);
   // bottom left split frame
   height = fSplitFrame->GetSecond()->GetFirst()->GetHeight();
   fSplitFrame->GetSecond()->GetFirst()->Resize(blwidth, height);
   // bottom center split frame
   height = fSplitFrame->GetSecond()->GetSecond()->GetFirst()->GetHeight();
   fSplitFrame->GetSecond()->GetSecond()->GetFirst()->Resize(bcwidth, height);
   // bottom right split frame
   height = fSplitFrame->GetSecond()->GetSecond()->GetSecond()->GetHeight();
   fSplitFrame->GetSecond()->GetSecond()->GetSecond()->Resize(brwidth, height);

   fSplitFrame->Layout();

   if (fIsEmbedded && gEve) {
      width = ((TGCompositeFrame *)gEve->GetBrowser()->GetTabBottom()->GetParent())->GetWidth();
      ((TGCompositeFrame *)gEve->GetBrowser()->GetTabBottom()->GetParent())->Resize(width, bottom_height);
      width = ((TGCompositeFrame *)gEve->GetBrowser()->GetTabRight()->GetParent())->GetWidth();
      ((TGCompositeFrame *)gEve->GetBrowser()->GetTabRight()->GetParent())->Resize(width, top_height);
   }
}

//______________________________________________________________________________
void SplitGLView::SaveConfig(const char *fname)
{

   Int_t bottom_height = 0;
   Int_t top_height = 0;
   TGSplitFrame *frm;
   TEnv *env = new TEnv(fname);

   if (fIsEmbedded && gEve) {
      env->SetValue("Eve.Width", (Int_t)gEve->GetBrowser()->GetWidth());
      env->SetValue("Eve.Height", (Int_t)gEve->GetBrowser()->GetHeight());
   }
   // get top (main) split frame
   frm = fSplitFrame->GetFirst();
   env->SetValue("MainView.Height", (Int_t)frm->GetHeight());
   // get bottom left split frame
   frm = fSplitFrame->GetSecond()->GetFirst();
   env->SetValue("Bottom.Left.Width", (Int_t)frm->GetWidth());
   // get bottom center split frame
   frm = fSplitFrame->GetSecond()->GetSecond()->GetFirst();
   env->SetValue("Bottom.Center.Width", (Int_t)frm->GetWidth());
   // get bottom right split frame
   frm = fSplitFrame->GetSecond()->GetSecond()->GetSecond();
   env->SetValue("Bottom.Right.Width", (Int_t)frm->GetWidth());
   if (fIsEmbedded && gEve) {
      top_height = (Int_t)((TGCompositeFrame *)gEve->GetBrowser()->GetTabRight()->GetParent())->GetHeight();
      env->SetValue("Right.Tab.Height", top_height);
      bottom_height = (Int_t)((TGCompositeFrame *)gEve->GetBrowser()->GetTabBottom()->GetParent())->GetHeight();
      env->SetValue("Bottom.Tab.Height", bottom_height);

      env->SetValue("Eve.Selection", gEve->GetSelection()->GetPickToSelect());
      env->SetValue("Eve.Highlight", gEve->GetHighlight()->GetPickToSelect());
   }

   env->SaveLevel(kEnvLocal);
#ifdef R__WIN32
   if (!gSystem->AccessPathName(Form("%s.new", fname))) {
      gSystem->Exec(Form("del %s", fname));
      gSystem->Rename(Form("%s.new", fname), fname);
   }
#endif
}

//______________________________________________________________________________
void SplitGLView::SwapToMainView(TGLViewerBase *viewer)
{
   // Swap frame embedded in a splitframe to the main view (slot method).

   TGCompositeFrame *parent = 0;
   if (!fSplitFrame->GetFirst()->GetFrame())
      return;
   if (viewer == 0) {
      TGPictureButton *src = (TGPictureButton*)gTQSender;
      parent = (TGCompositeFrame *)src->GetParent();
      while (parent && !parent->InheritsFrom("TGSplitFrame")) {
         parent = (TGCompositeFrame *)parent->GetParent();
      }
   }
   else {
      TGCompositeFrame *src = ((TGLEmbeddedViewer *)viewer)->GetFrame();
      if (!src) return;
      TGLOverlayButton *but = (TGLOverlayButton *)((TQObject *)gTQSender);
      but->ResetState();
      parent = (TGCompositeFrame *)src->GetParent();
   }
   if (parent && parent->InheritsFrom("TGSplitFrame"))
      ((TGSplitFrame *)parent)->SwitchToMain();
}

//______________________________________________________________________________
void SplitGLView::UnDock(TGLViewerBase *viewer)
{
   // Undock frame embedded in a splitframe (slot method).

   TGCompositeFrame *src = ((TGLEmbeddedViewer *)viewer)->GetFrame();
   if (!src) return;
   TGLOverlayButton *but = (TGLOverlayButton *)((TQObject *)gTQSender);
   but->ResetState();
   TGCompositeFrame *parent = (TGCompositeFrame *)src->GetParent();
   if (parent && parent->InheritsFrom("TGSplitFrame"))
      ((TGSplitFrame *)parent)->ExtractFrame();
}

//______________________________________________________________________________
void SplitGLView::UpdateSummary()
{
   // Update summary of current event.

   TEveElement::List_i i;
   TEveElement::List_i j;
   Int_t k;
   TEveElement *el;
   HtmlObjTable *table;
   TEveEventManager *mgr = gEve ? gEve->GetCurrentEvent() : 0;
   if (mgr) {
      fgHtmlSummary->Clear("D");
      for (i=mgr->BeginChildren(); i!=mgr->EndChildren(); ++i) {
         el = ((TEveElement*)(*i));
         if (el->IsA() == TEvePointSet::Class()) {
            TEvePointSet *ps = (TEvePointSet *)el;
            TString ename  = ps->GetElementName();
            TString etitle = ps->GetElementTitle();
            if (ename.First('\'') != kNPOS)
               ename.Remove(ename.First('\''));
            etitle.Remove(0, 2);
            Int_t nel = atoi(etitle.Data());
            table = fgHtmlSummary->AddTable(ename, 0, nel);
         }
         else if (el->IsA() == TEveTrackList::Class()) {
            TEveTrackList *tracks = (TEveTrackList *)el;
            TString ename  = tracks->GetElementName();
            if (ename.First('\'') != kNPOS)
               ename.Remove(ename.First('\''));
            table = fgHtmlSummary->AddTable(ename.Data(), 5,
                     tracks->NumChildren(), kTRUE, "first");
            table->SetLabel(0, "Momentum");
            table->SetLabel(1, "P_t");
            table->SetLabel(2, "Phi");
            table->SetLabel(3, "Theta");
            table->SetLabel(4, "Eta");
            k=0;
            for (j=tracks->BeginChildren(); j!=tracks->EndChildren(); ++j) {
               Float_t p     = ((TEveTrack*)(*j))->GetMomentum().Mag();
               table->SetValue(0, k, p);
               Float_t pt    = ((TEveTrack*)(*j))->GetMomentum().Perp();
               table->SetValue(1, k, pt);
               Float_t phi   = ((TEveTrack*)(*j))->GetMomentum().Phi();
               table->SetValue(2, k, phi);
               Float_t theta = ((TEveTrack*)(*j))->GetMomentum().Theta();
               table->SetValue(3, k, theta);
               Float_t eta   = ((TEveTrack*)(*j))->GetMomentum().Eta();
               table->SetValue(4, k, eta);
               ++k;
            }
         }
      }
      fgHtmlSummary->Build();
      fgHtml->Clear();
      fgHtml->ParseText((char*)fgHtmlSummary->Html().Data());
      fgHtml->Layout();
   }
}


#ifdef __CINT__
void SplitGLView()
{
   printf("This script is used via ACLiC by the macro \"alice_esd_split.C\"\n");
   printf("To see it in action, just run \".x alice_esd_split.C\"\n");
   return;
}
#endif


