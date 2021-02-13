/*
 * $Header$
 * $Log$
 */

#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>

#include <iostream>
#include <TSystem.h>
#include <TGTab.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGWindow.h>
#include <TGButton.h>
#include <TGMsgBox.h>
#include <TVirtualX.h>
#include <strlcpy.h>
#include <snprintf.h>

#include "XSVarious.h"
#include "XSElementDlg.h"
#include "XSReactionDlg.h"

#include "NdbMTDir.h"

//ClassImp(XSReactionDlg);

// Options for filling containers
enum   DirOptions {
      DIROnlyFiles,
      DIROnlyDirectories,
      DIRBoth
   };

// Options for creating the path
enum   PathOptions {
      PATHIsotope,
      PATHProjectile,
      PATHDatabase,
      PATHFile
   };

enum   XSReactionMessages {
      REAC_ELEMENT_STEP,
      REAC_TABLE,
      REAC_OK,
      REAC_EXEC,
      REAC_CLOSE,
      REAC_RESET,
      REAC_PTABLE
   };

enum   ComboIds {
      ISOTOPE_COMBO,
      PROJECTILE_COMBO,
      TEMPERATURE_COMBO,
      DATABASE_COMBO,
      LINE_WIDTH_COMBO,
      LINE_COLOR_COMBO,
      MARKER_STYLE_COMBO,
      MARKER_COLOR_COMBO,
      MARKER_SIZE_COMBO,
      ERRORBAR_COLOR_COMBO,
      REACTION_LISTBOX,
   };

static   Int_t   LastWinX = -1;
static   Int_t   LastWinY = -1;

/* ----- XSReactionDlg ----- */
XSReactionDlg::XSReactionDlg( const TGWindow *p,
      const TGWindow *main, UInt_t initZ, UInt_t w, UInt_t h)
   : TGTransientFrame(p,main,w,h)
{
   // Remember the main window
   mainWindow = main;

   /* ---------- Prepare the Layout Hints ------- */
   lHFixed = new TGLayoutHints(kLHintsTop | kLHintsLeft,
               3, 3, 2, 2);
   lHExpX = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
               3, 3, 2, 2);
   lHExpY = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandY,
               3, 3, 2, 2);
   lHExpXY = new TGLayoutHints(kLHintsTop | kLHintsLeft |
            kLHintsExpandX | kLHintsExpandY,
               3, 3, 2, 2);
   lHBot = new TGLayoutHints(kLHintsBottom | kLHintsRight,
               3, 3, 2, 2);

   lHExpXCen = new TGLayoutHints( kLHintsTop | kLHintsLeft |
            kLHintsExpandX | kLHintsCenterY,
            5, 5, 2, 2);
   lHFixedCen = new TGLayoutHints( kLHintsTop | kLHintsLeft |
            kLHintsCenterY,
            5, 5, 2, 2);

   /* --------- Create the Dialog -------------- */
   buttonFrame = new TGHorizontalFrame(this, 240, 20, kFixedWidth);
   okButton = new TGTextButton(buttonFrame, "&Ok", REAC_OK);
   okButton->Associate(this);
   execButton = new TGTextButton(buttonFrame, "&Execute", REAC_EXEC);
   execButton->Associate(this);
   resetButton = new TGTextButton(buttonFrame, "&Reset", REAC_RESET);
   resetButton->Associate(this);
   closeButton = new TGTextButton(buttonFrame, "&Close", REAC_CLOSE);
   closeButton->Associate(this);

   buttonFrame->AddFrame(okButton, lHExpX);
   buttonFrame->AddFrame(execButton, lHExpX);
   buttonFrame->AddFrame(resetButton, lHExpX);
   buttonFrame->AddFrame(closeButton, lHExpX);
   AddFrame(buttonFrame, lHBot);

   /* ------------ Material --------------- */
   materialGroup = new TGGroupFrame(this,"Material",kVerticalFrame);

   /* ----- First Sub-Frame ---- */
   frm1 = new TGHorizontalFrame(
               materialGroup, w, 36,
               kChildFrame|kFitWidth );
   frm1->ChangeBackground(0xFF00);

   elementLbl = new TGLabel(frm1,"Element:");
   elementLbl->SetTextJustify(kTextRight | kTextCenterY);
   frm1->AddFrame(elementLbl,lHFixedCen);

   elementBuf = new TGTextBuffer(20);
   elementText = new TGTextEntry(frm1,elementBuf);
   elementText->Resize(50, elementText->GetDefaultHeight());
   elementText->Associate(this);
   frm1->AddFrame(elementText,lHFixedCen);

   elementStep = new XSStepButton(frm1,REAC_ELEMENT_STEP);
   elementStep->Associate(this);
   frm1->AddFrame(elementStep,lHFixedCen);

   nameLbl = new TGLabel(frm1,"X",blueBoldGC);
   nameLbl->SetTextJustify(kTextLeft | kTextCenterY);
   frm1->AddFrame(nameLbl,lHExpXCen);

   ptableButton = new TGPictureButton(frm1,
            fClient->GetPicture(PTBL_ICON),
            REAC_PTABLE);
   ptableButton->SetToolTipText(
      "Choose the Element from the Periodic Table");
   ptableButton->Associate(this);
   frm1->AddFrame(ptableButton,lHFixedCen);

   materialGroup->AddFrame(frm1,lHExpX);

   /* ----- Second Sub-Frame ---- */
   frm2 = new TGHorizontalFrame(materialGroup,w,20,
            kChildFrame|kFitWidth);

   Vfrm1 = new TGCompositeFrame(frm2, w, 70,
               kChildFrame|kFitHeight|kFitWidth);
   Vfrm1->SetLayoutManager(new TGMatrixLayout(Vfrm1,0,4));

   chargeLbl = new TGLabel(Vfrm1,"Charge (Z):");
   chargeLbl->SetTextJustify(kTextRight | kTextCenterY);
   Vfrm1->AddFrame(chargeLbl);

   zLbl = new TGLabel(Vfrm1,"ZZZZZZZZ",blueBoldGC);
   zLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm1->AddFrame(zLbl);

   massLbl = new TGLabel(Vfrm1,"Atomic Mass:");
   massLbl->SetTextJustify(kTextRight | kTextCenterY);
   Vfrm1->AddFrame(massLbl);

   massValLbl = new TGLabel(Vfrm1,"ZZZZZZZZ",blueBoldGC);
   massValLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm1->AddFrame(massValLbl);

   densityLbl = new TGLabel(Vfrm1,"Density:");
   densityLbl->SetTextJustify(kTextRight | kTextCenterY);
   Vfrm1->AddFrame(densityLbl);

   densityValLbl = new TGLabel(Vfrm1,"ZZZZZZZZ",blueBoldGC);
   densityValLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm1->AddFrame(densityValLbl);

   oxidationLbl = new TGLabel(Vfrm1,"Oxidation:");
   oxidationLbl->SetTextJustify(kTextRight | kTextCenterY);
   Vfrm1->AddFrame(oxidationLbl);

   oxidationValLbl = new TGLabel(Vfrm1,"ZZZZZZZZ",blueBoldGC);
   oxidationValLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm1->AddFrame(oxidationValLbl);

   meltingPtLbl = new TGLabel(Vfrm1,"Melting Pt (C):");
   meltingPtLbl->SetTextJustify(kTextRight | kTextCenterY);
   Vfrm1->AddFrame(meltingPtLbl);

   meltingValLbl = new TGLabel(Vfrm1,"ZZZZZZZZ",blueBoldGC);
   meltingValLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm1->AddFrame(meltingValLbl);

   boilingPtLbl = new TGLabel(Vfrm1,"Boiling Pt (C):");
   boilingPtLbl->SetTextJustify(kTextRight | kTextCenterY);
   Vfrm1->AddFrame(boilingPtLbl);

   boilingValLbl = new TGLabel(Vfrm1,"ZZZZZZZZ",blueBoldGC);
   boilingValLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm1->AddFrame(boilingValLbl);

   frm2->AddFrame(Vfrm1,lHExpX);
   materialGroup->AddFrame(frm2,lHExpX);

   /* ----- Fourth Sub-Frame ---- */
   frm3 = new TGHorizontalFrame(materialGroup,w,20,
            kChildFrame|kFitWidth);

   isotopeLbl = new TGLabel(frm3,"Isotope (A):");
   isotopeLbl->SetTextJustify(kTextRight | kTextCenterY);
   frm3->AddFrame(isotopeLbl,lHFixedCen);

   isotopeCombo = new TGComboBox(frm3,ISOTOPE_COMBO);
   isotopeCombo->Resize(80,20);
   isotopeCombo->Associate(this);
   frm3->AddFrame(isotopeCombo,lHFixedCen);

   isotopeInfoLbl = new TGLabel(frm3,"Isotope Info:");
   isotopeInfoLbl->SetTextJustify(kTextRight | kTextCenterY);
   frm3->AddFrame(isotopeInfoLbl,lHFixedCen);

   isotopeInfoValLbl = new TGLabel(frm3,"MMMM",blueBoldGC);
   isotopeInfoValLbl->SetTextJustify(kTextLeft | kTextCenterY);
   frm3->AddFrame(isotopeInfoValLbl,lHExpXCen);

   materialGroup->AddFrame(frm3,lHExpX);

   AddFrame(materialGroup,lHExpX);

   /* --------------- Reaction ------------------ */
   reactionGroup = new TGGroupFrame(this,"Reaction",kVerticalFrame);

   frm4 = new TGHorizontalFrame(reactionGroup,w,60,
               kChildFrame|kFitWidth);

   Vfrm2 = new TGCompositeFrame(frm4, 160, 70,
               kChildFrame|kFitHeight|kFitWidth);

   Vfrm2->SetLayoutManager(new TGMatrixLayout(Vfrm2,0,2));

   projectileLbl = new TGLabel(Vfrm2,"Projectile:");
   projectileLbl->SetTextJustify(kTextRight | kTextCenterY);
   Vfrm2->AddFrame(projectileLbl);

   projectileCombo = new TGComboBox(Vfrm2,PROJECTILE_COMBO);
   projectileCombo->Resize(90,20);
   projectileCombo->Associate(this);
   Vfrm2->AddFrame(projectileCombo);

   temperatureLbl = new TGLabel(Vfrm2,"Temperature:");
   temperatureLbl->SetTextJustify(kTextRight | kTextCenterY);
   Vfrm2->AddFrame(temperatureLbl);

   temperatureCombo = new TGComboBox(Vfrm2,TEMPERATURE_COMBO);
   temperatureCombo->Resize(90,20);
   temperatureCombo->Associate(this);
   Vfrm2->AddFrame(temperatureCombo);

   databaseLbl = new TGLabel(Vfrm2,"Database:");
   databaseLbl->SetTextJustify(kTextRight | kTextCenterY);
   Vfrm2->AddFrame(databaseLbl);

   databaseCombo = new TGComboBox(Vfrm2,DATABASE_COMBO);
   databaseCombo->Resize(90,20);
   databaseCombo->Associate(this);
   Vfrm2->AddFrame(databaseCombo);

   frm4->AddFrame(Vfrm2,lHExpXCen);

   reactionLbl = new TGLabel(frm4,"Reaction:");
   reactionLbl->SetTextJustify(kTextLeft | kTextCenterY);
   frm4->AddFrame(reactionLbl,lHFixedCen);

   reactionList = new TGListBox(frm4,REACTION_LISTBOX);
   reactionList->Resize(70,80);
   //reactionList->SetMultipleSelections(kTRUE);
   reactionList->Associate(this);
   frm4->AddFrame(reactionList,lHExpXCen);

   reactionGroup->AddFrame(frm4,lHExpX);

   // --- Second sub frame ---
   frm5 = new TGHorizontalFrame(reactionGroup,w,20,
               kChildFrame|kFitWidth);
   reactionInfoLbl = new TGLabel(frm5,"Reaction Info:");
   reactionInfoLbl->SetTextJustify(kTextRight | kTextCenterY);
   frm5->AddFrame(reactionInfoLbl,lHFixedCen);

   reactionInfoValLbl = new TGLabel(frm5,"-",blueBoldGC);
   reactionInfoValLbl->SetTextJustify(kTextLeft | kTextCenterY);
   frm5->AddFrame(reactionInfoValLbl,lHExpXCen);

   reactionGroup->AddFrame(frm5,lHExpX);

   AddFrame(reactionGroup,lHExpX);

   /* --------------- Options ------------------ */
   optionGroup = new TGGroupFrame(this,"Options",kVerticalFrame);

   Vfrm3 = new TGCompositeFrame(optionGroup, w, 70,
            kChildFrame|kFitHeight|kFitWidth);

   Vfrm3->SetLayoutManager(new TGMatrixLayout(Vfrm3,0,4));

   lineWidthLbl = new TGLabel(Vfrm3,"Line Width:");
   lineWidthLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm3->AddFrame(lineWidthLbl);

   lineWidthCombo = new TGComboBox(Vfrm3,LINE_WIDTH_COMBO);
   lineWidthCombo->Resize(100,20);
   lineWidthCombo->Associate(this);
   Vfrm3->AddFrame(lineWidthCombo);

   lineColorLbl = new TGLabel(Vfrm3,"    Line Color:");
   lineColorLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm3->AddFrame(lineColorLbl);

   lineColorCombo = new TGComboBox(Vfrm3,LINE_COLOR_COMBO);
   lineColorCombo->Resize(100,20);
   lineColorCombo->Associate(this);
   Vfrm3->AddFrame(lineColorCombo);

   markerStyleLbl = new TGLabel(Vfrm3,"Marker Style:");
   markerStyleLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm3->AddFrame(markerStyleLbl);

   markerStyleCombo = new TGComboBox(Vfrm3,MARKER_STYLE_COMBO);
   markerStyleCombo->Resize(100,20);
   markerStyleCombo->Associate(this);
   Vfrm3->AddFrame(markerStyleCombo);

   markerColorLbl = new TGLabel(Vfrm3,"    Marker Color:");
   markerColorLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm3->AddFrame(markerColorLbl);

   markerColorCombo = new TGComboBox(Vfrm3, MARKER_COLOR_COMBO);
   markerColorCombo->Resize(100,20);
   markerColorCombo->Associate(this);
   Vfrm3->AddFrame(markerColorCombo);

   markerSizeLbl = new TGLabel(Vfrm3,"Marker Size:");
   markerSizeLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm3->AddFrame(markerSizeLbl);

   markerSizeCombo = new TGComboBox(Vfrm3, MARKER_SIZE_COMBO);
   markerSizeCombo->Resize(100,20);
   markerSizeCombo->Associate(this);
   Vfrm3->AddFrame(markerSizeCombo);

   errorbarColorLbl = new TGLabel(Vfrm3,"    ErrorBar Color:");
   errorbarColorLbl->SetTextJustify(kTextLeft | kTextCenterY);
   Vfrm3->AddFrame(errorbarColorLbl);

   errorbarColorCombo = new TGComboBox(Vfrm3,ERRORBAR_COLOR_COMBO);
   errorbarColorCombo->Resize(100,20);
   errorbarColorCombo->Associate(this);
   Vfrm3->AddFrame(errorbarColorCombo);

   optionGroup->AddFrame(Vfrm3,lHExpX);
   AddFrame(optionGroup,lHExpX);

   /* -------------- Information ------------------ */
   infoGroup = new TGGroupFrame(this,"Information",kVerticalFrame);

   infoView = new TGTextView(infoGroup, 200, 100,
            kChildFrame | kSunkenFrame);
   infoGroup->AddFrame(infoView,lHExpXY);

   AddFrame(infoGroup, lHExpXY);

   /* --------- Initialize layout algorithm with Resize() -------- */
   MapSubwindows();
   Resize(GetDefaultSize());

   /* --------- Initialise with Element Z ---------- */
   SetElement(initZ);
   InitCombos();

   /* --------- Set Windows Position --------- */
   int   ax, ay;

   if (LastWinX == -1) {   // Go to the middle of parent window
      Window_t wdum;
      gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
         (((TGFrame *) main)->GetWidth() - fWidth) >> 1,
         (((TGFrame *) main)->GetHeight() - fHeight) >> 1,
         ax, ay, wdum);
      ax = MAX(ax,0);
      ay = MAX(ay,0);
   } else {
      ax = LastWinX;
      ay = LastWinY;
   }
   Move(ax,ay);
   SetWMPosition(ax,ay);

   /* --------- Map everything -------- */
   SetWindowName("Select Element/Reaction");
   MapWindow();

} // XSReactionDlg

/* ----- ~XSReactionDlg ----- */
XSReactionDlg::~XSReactionDlg()
{
   // --- Delete Layout Hints ---
   delete   lHFixed;
   delete   lHExpX;
   delete   lHExpY;
   delete   lHExpXY;
   delete   lHFixedCen;
   delete   lHExpXCen;

   // --- Delete Material Group ----
   delete   elementLbl;
   delete   elementText;
   delete   elementStep;
   delete   ptableButton;
   delete   chargeLbl;
   delete   zLbl;
   delete   nameLbl;
   delete   massLbl;
   delete   massValLbl;
   delete   densityLbl;
   delete   densityValLbl;
   delete   meltingPtLbl;
   delete   meltingValLbl;
   delete   boilingPtLbl;
   delete   boilingValLbl;
   delete   oxidationLbl;
   delete   oxidationValLbl;
   delete   isotopeInfoLbl;
   delete   isotopeInfoValLbl;
   delete   isotopeLbl;
   delete   isotopeCombo;
   delete   materialGroup;

   // --- Delete Reaction group ---
   delete   projectileLbl;
   delete   projectileCombo;
   delete   temperatureLbl;
   delete   temperatureCombo;
   delete   databaseLbl;
   delete   databaseCombo;
   delete   reactionLbl;
   delete   reactionList;
   delete   reactionInfoLbl;
   delete   reactionInfoValLbl;
   delete   reactionGroup;

   // --- Option Group ----
   delete   lineWidthLbl;
   delete   lineColorLbl;
   delete   markerStyleLbl;
   delete   markerColorLbl;
   delete   markerSizeLbl;
   delete   errorbarColorLbl;

   delete   lineWidthCombo;
   delete   lineColorCombo;
   delete   markerStyleCombo;
   delete   markerColorCombo;
   delete   markerSizeCombo;
   delete   errorbarColorCombo;

   delete   optionGroup;

   // --- Info Group ----
   delete   infoView;
   delete   infoGroup;

   // --- Frames ---
   delete   frm1;
   delete   frm2;
   delete   frm3;
   delete   frm4;
   delete   frm5;

   delete   Vfrm1;
   delete   Vfrm2;
   delete   Vfrm3;

   // --- Buttons ---
   delete   okButton;
   delete   execButton;
   delete   resetButton;
   delete   closeButton;
   delete   buttonFrame;
} // ~XSReactionDlg

/* ----- InitColorCombo ----- */
void
XSReactionDlg::InitColorCombo(TGComboBox *cb)
{
   // Normally this should be filled with color entries!!!

   cb->AddEntry("Black",0);
   cb->AddEntry("Red",1);
   cb->AddEntry("Green",2);
   cb->AddEntry("Blue",3);
   cb->AddEntry("Yellow",4);
   cb->AddEntry("Magenta",5);
   cb->AddEntry("Cyan",6);
   // .......
   cb->Select(0);      // <<<< Should find the correct color...
} // InitColorCombos

/* ----- InitCombos ----- */
void
XSReactionDlg::InitCombos()
{
   InitColorCombo(lineColorCombo);
   InitColorCombo(markerColorCombo);
   InitColorCombo(errorbarColorCombo);

   /* ---- Line Widths ----- */
   lineWidthCombo->AddEntry("None",0);
   lineWidthCombo->AddEntry("1",1);
   lineWidthCombo->AddEntry("2",2);
   lineWidthCombo->AddEntry("3",3);
   lineWidthCombo->AddEntry("4",4);
   lineWidthCombo->AddEntry("5",5);
   lineWidthCombo->AddEntry("6",6);
   lineWidthCombo->AddEntry("7",7);
   lineWidthCombo->AddEntry("8",8);
   lineWidthCombo->AddEntry("9",9);
   lineWidthCombo->AddEntry("10",10);
   lineWidthCombo->Select(1);
   // --- Select the correct one ---

   /* ---- Marker Style ----- */
   markerStyleCombo->AddEntry("None",0);
   markerStyleCombo->AddEntry("Bullet",1);
   markerStyleCombo->AddEntry("Triangle Up",2);
   markerStyleCombo->AddEntry("Triangle Down",3);
   markerStyleCombo->AddEntry("Square",4);
   markerStyleCombo->AddEntry("Diamond",5);
   markerStyleCombo->AddEntry("Star",6);
   markerStyleCombo->AddEntry("Empty Bullet",7);
   markerStyleCombo->AddEntry("Empty Tri-Up",8);
   markerStyleCombo->AddEntry("Empty Tri-Down",9);
   markerStyleCombo->AddEntry("Empty Square",10);
   markerStyleCombo->Select(0);
   // --- Select the correct one ---

   /* ---- Marker Size ----- */
   markerSizeCombo->AddEntry("None",0);
   markerSizeCombo->AddEntry("0.1",1);
   markerSizeCombo->AddEntry("0.2",2);
   markerSizeCombo->AddEntry("0.3",3);
   markerSizeCombo->AddEntry("0.4",4);
   markerSizeCombo->AddEntry("0.5",5);
   markerSizeCombo->AddEntry("0.6",6);
   markerSizeCombo->AddEntry("0.7",7);
   markerSizeCombo->AddEntry("0.8",8);
   markerSizeCombo->AddEntry("0.9",9);
   markerSizeCombo->AddEntry("1.0",10);
   markerSizeCombo->Select(0);

   // --- Select the correct one ---
} // InitCombos

/* ----- GetString ------- */
const char *
XSReactionDlg::GetString( int box )
{
   const TGLBEntry   *entry;
   switch  (box) {
      case ISOTOPE_COMBO:
         entry = isotopeCombo->GetSelectedEntry();
         break;

      case PROJECTILE_COMBO:
         entry = projectileCombo->GetSelectedEntry();
         break;

      case TEMPERATURE_COMBO:
         entry = temperatureCombo->GetSelectedEntry();
         break;

      case DATABASE_COMBO:
         entry = databaseCombo->GetSelectedEntry();
         break;

      case LINE_WIDTH_COMBO:
         entry = lineWidthCombo->GetSelectedEntry();
         break;

      case LINE_COLOR_COMBO:
         entry = lineColorCombo->GetSelectedEntry();
         break;

      case MARKER_STYLE_COMBO:
         entry = markerStyleCombo->GetSelectedEntry();
         break;

      case MARKER_COLOR_COMBO:
         entry = markerColorCombo->GetSelectedEntry();
         break;

      case MARKER_SIZE_COMBO:
         entry = markerSizeCombo->GetSelectedEntry();
         break;

      case ERRORBAR_COLOR_COMBO:
         entry = errorbarColorCombo->GetSelectedEntry();
         break;

      case REACTION_LISTBOX:
         entry = reactionList->GetSelectedEntry();
         break;

      default:
         entry = NULL;
   }
   if (entry==NULL) return NULL;

   return ((TGTextLBEntry*)entry)->GetText()->GetString();
} // GetString

/* ----- CreatePath ------ */
char*
XSReactionDlg::CreatePath( int option )
{
   static   char   path[256];

   strlcpy(path, DBDIR,256);   // Initialise directory
   strlcat(path, XSelements->Mnemonic(Z),256);
   strlcat(path, PATHSEP,256);

   if (option == PATHIsotope) return path;

   /* --- add the selected isotope --- */
   strlcat(path, GetString(ISOTOPE_COMBO),256);
   strlcat(path, PATHSEP,256);

   if (option == PATHProjectile) return path;

   /* --- add the selected projectile --- */
   strlcat(path, GetString(PROJECTILE_COMBO),256);
   strlcat(path, PATHSEP,256);

   if (option == PATHDatabase) return path;

   /* --- finally add the file --- */
   strlcat(path, GetString(DATABASE_COMBO),256);

   return path;
} // CreatePath

/* ----- UpdateContainer ----- */
int
XSReactionDlg::UpdateContainer( TGListBox *lb, char *path, int option)
{
   const char *entry;
   FileStat_t st;

   // --- First Remove everything ---
   lb->RemoveEntries(0,1000);

   // Scan directory to update the combo box
   void *dirp = gSystem->OpenDirectory(path);
   if (dirp==NULL) {
      lb->AddEntry("-",0);
      return kTRUE;   // Ooops not found
   }

   int   i=0;
   while ((entry = gSystem->GetDirEntry(dirp))) {
      // Skip the . and .* directory entries
      if (entry[0] == '.') continue;

      char   fn[256];
      strlcpy(fn,path,256);
      strlcat(fn,entry,256);

      gSystem->GetPathInfo(fn, st);
      if (((option == DIROnlyFiles) && !R_ISDIR(st.fMode)) ||
          ((option == DIROnlyDirectories) && R_ISDIR(st.fMode)) ||
          (option == DIRBoth)) {
         lb->AddEntry(entry,i++);
      }
   }
   gSystem->FreeDirectory(dirp);
   if (i==0) {
      lb->AddEntry("None",0);
      return kTRUE;   // Ooops not found
   }
   return kFALSE;
} // UpdateContainer

/* ----- UpdateIsotopes ------ */
void
XSReactionDlg::UpdateIsotopes()
{
   // --- Update the isotopes ---
   UpdateContainer(
      isotopeCombo->GetListBox(),
      CreatePath(PATHIsotope),
      DIROnlyDirectories);
   isotopeCombo->Select(0);
} // UpdateIsotopes

/* ----- UpdateProjectile ------ */
void
XSReactionDlg::UpdateProjectile()
{
   // --- Update the projectiles ---
   UpdateContainer(
      projectileCombo->GetListBox(),
      CreatePath(PATHProjectile),
      DIROnlyDirectories);
   projectileCombo->Select(0);
} // UpdateProjectile

/* ----- UpdateDatabase ------ */
void
XSReactionDlg::UpdateDatabase()
{
   // --- Update the databases ---
   UpdateContainer(
      databaseCombo->GetListBox(),
      CreatePath(PATHDatabase),
      DIRBoth);
   databaseCombo->Select(0);
} // UpdateDatabase

/* ----- UpdateReactions ------ */
void
XSReactionDlg::UpdateReactions()
{
   NdbMTDir   elemDir;

   reactionList->RemoveEntries(0,1000);
   infoView->Clear();

   if (elemDir.LoadENDF(CreatePath(PATHFile)))
      return;

   int id = 0;
   for (int i=0; i<elemDir.Sections(); i++) {
      if (elemDir.DIRMF(i) == 3)
         reactionList->AddEntry(
            XSReactionDesc->GetShort(elemDir.DIRMT(i)),
            id++);
   }
        reactionList->MapSubwindows();
   reactionList->Select(0);

   // ------------------------------------------
   // Prepare a string with information
   TString   info;

   info.Append("Symbol name:\t\t");
   info.Append(elemDir.SymbolName());

   info.Append("\nLaboratory:\t\t");
   info.Append(elemDir.Laboratory());

   info.Append("\nEvaluation Date:\t");
   info.Append(elemDir.EvaluationDate());

   info.Append("\nAuthor(s):\t\t");
   info.Append(elemDir.Author());

   info.Append("\nReference:\t\t");
   info.Append(elemDir.Reference());

   info.Append("\nDistribution Date:\t");
   info.Append(elemDir.DistributionDate());
   info.Append("\nLast Revision Date:\t");
   info.Append(elemDir.LastRevisionDate());
   info.Append("\nMaster Entry Date:\t");
   info.Append(elemDir.MasterEntryDate());
   info.Append("\n\n");
   info.Append(elemDir.GetInfo());

   infoView->LoadBuffer(info.Data());
} // UpdateReactions

/* ----- UpdateCurIsotope ----- */
void
XSReactionDlg::UpdateCurIsotope()
{
   XSElement   *elem = XSelements->Elem(Z);

   const char *str = GetString(ISOTOPE_COMBO);
   if (str)
      isotopeInfoValLbl->SetText(
         new TGString(elem->IsotopeInfo(str)));
   else
      isotopeInfoValLbl->SetText(new TGString("-"));

   // --- Update combos depending on isotope ---
   UpdateProjectile();
   UpdateDatabase();
   UpdateReactions();
} // UpdateCurIsotope

/* ----- SetElement ----- */
void
XSReactionDlg::SetElement(UInt_t aZ)
{
   char   str[5];

   Z = aZ;
   if (Z<1)
      Z = 1;
   else
   if (Z>XSelements->GetSize())
      Z = XSelements->GetSize();


   XSElement   *elem = XSelements->Elem(Z);

   elementBuf->Clear();
   elementBuf->AddText(0,elem->Symbol());
   fClient->NeedRedraw(elementText);
   nameLbl->SetText(new TGString(elem->Name()));

   snprintf(str,5,"%d",Z);
   zLbl->SetText(new TGString(str));

   // --- Update several values for element ---
   massValLbl->SetText(new TGString(elem->AtomicWeight()));
   densityValLbl->SetText(new TGString(elem->Density()));
   densityValLbl->SetText(new TGString(elem->Density()));
   meltingValLbl->SetText(new TGString(elem->MeltingPt()));
   boilingValLbl->SetText(new TGString(elem->BoilingPt()));
   oxidationValLbl->SetText(new TGString(elem->Oxidation()));

   // --------------------------------
   UpdateIsotopes();
   UpdateCurIsotope();

   // --------------------------------
   // ---Update the rest of combos ---
   // --------------------------------
   temperatureCombo->RemoveEntries(0,1000);
   temperatureCombo->AddEntry("0",0);
   temperatureCombo->AddEntry("293",1);
   temperatureCombo->AddEntry("900",2);
   temperatureCombo->AddEntry("1000",3);
   temperatureCombo->Select(1);


   // Redraw everything
   fClient->NeedRedraw(this);
   Layout();
} // SetElement

/* ----- CloseWindow ----- */
void
XSReactionDlg::CloseWindow()
{
   // --- Remember old position ---
   Window_t wdum;
   gVirtualX->TranslateCoordinates(GetId(), GetParent()->GetId(),
      0, 0, LastWinX, LastWinY, wdum);

   delete this;
} // CloseWindow

/* ----- elementEntryChanged ---- */
void
XSReactionDlg::ElementEntryChanged()
{
   const   char   *name;
   UInt_t   newZ;

   /* --- get contents --- */
   name = elementBuf->GetString();

   /* --- check if it is a number --- */
   newZ = atoi(name);
   if ((newZ>0) && (newZ<=XSelements->GetSize())) {
      SetElement(newZ);
      return;
   }

   /* Search for the specific element */
   newZ = XSelements->Find(name);

   if (newZ == 0) {
      // --- Error Entry ---
      // issue message and put the old one
      char   msg[256];
      int   retval;
      snprintf(msg,256,"Element %s not found!", elementBuf->GetString());
      new TGMsgBox(fClient->GetRoot(),this,"Not found",msg,
            kMBIconAsterisk,kMBOk,&retval);
      SetElement(Z);
   } else
      SetElement(newZ);
} // ElementEntryChanged

/* ----- UpdateGraph ---- */
void
XSReactionDlg::UpdateGraph( NdbMTReactionXS *xs )
{
   char title[256];

   //!!! --- Only if new ---
   //!!! --- Else update the old one ---
   XSGraph *gr = new XSGraph(xs);

   // Prepare the title
   snprintf(title,256,"%s-%s %s",
         XSelements->Mnemonic(Z),
         GetString(ISOTOPE_COMBO),
         xs->Description().Data()
      );

   gr->GetGraph()->SetTitle(title);
   gr->GetGraph()->SetFillColor(19);
   gr->GetGraph()->SetLineColor(2);
   gr->GetGraph()->SetLineWidth(2);
   gr->GetGraph()->Draw("ALW");

//   gr->GetGraph()->GetHistogram()->SetXTitle(TString("Energy (eV)"));
//   gr->GetGraph()->GetHistogram()->SetYTitle(TString("Cross Section (barn)"));

   Add2GraphList(gr);
   fClient->NeedRedraw(canvasWindow);
   fClient->NeedRedraw(const_cast<TGWindow*>(mainWindow));
   canvas->Update();
} // UpdateGraph

/* ----- ExecCommand ---- */
Bool_t
XSReactionDlg::ExecCommand()
{
   // --- Load the data base ---
   NdbMTDir   el;
   if (el.LoadENDF(CreatePath(PATHFile)))
      return kTRUE;

   // --- Get the selected reaction ---
   const char *reacstr = GetString(REACTION_LISTBOX);

   // --- Search for the reaction ---
   for (int i=0; i<el.Sections(); i++) {
      if (el.DIRMF(i) == 3 && XSReactionDesc->GetShort(el.DIRMT(i))) {
         if (reacstr && !strcmp(reacstr,
              XSReactionDesc->GetShort(el.DIRMT(i)))) {
            NdbMTReactionXS   xs(el.DIRMT(i), reacstr);
            xs.LoadENDF(CreatePath(PATHFile));

            UpdateGraph(&xs);
         }
      }
   }
   return kFALSE;
} // ExecCommand

/* ----- ProcessButton ----- */
Bool_t
XSReactionDlg::ProcessButton(Long_t param1, Long_t param2)
{
   UInt_t   newZ;
   switch (param1) {
      case REAC_ELEMENT_STEP:
         if (param2==XSSTEPBUTTON_UP) {
            newZ = Z+1;
            if (newZ > XSelements->GetSize()) Z = 1;
         } else {
            newZ = Z-1;
            if (newZ<1) Z = XSelements->GetSize();
         }
         SetElement(newZ);
         break;

      case REAC_OK:
         // Execute the command and close the window
         if (!ExecCommand())
            CloseWindow();
         break;

      case REAC_EXEC:
         // Execute the command, but don't close the window
         ExecCommand();
         break;

      case REAC_RESET:
         SetElement(1);
         break;

      case REAC_CLOSE:
         CloseWindow();
         break;

      case REAC_PTABLE:
         newZ = Z;
         new XSElementDlg(fClient->GetRoot(), this, &newZ);
         if (newZ>0)
            SetElement(newZ);
         break;

      default:
         break;
   }
   return kTRUE;
} // ProcessButton

/* ----- ProcessCombo ------ */
Bool_t
XSReactionDlg::ProcessCombo( Long_t param1, Long_t param2 )
{
   printf("ComboMessage %ld %ld\n",param1,param2);

   switch (param1) {
      case ISOTOPE_COMBO:
         UpdateCurIsotope();
         break;

      case PROJECTILE_COMBO:
         UpdateDatabase();
         UpdateReactions();
         break;

      case DATABASE_COMBO:
         UpdateReactions();
         break;

      case REACTION_LISTBOX:
//////         reactionInfoValLbl->SetText(
//////            new TGString(elem->IsotopeInfo(str)));
         break;

      default:
         // Ooops an error in the code
         break;
   }

   return kTRUE;
} // ProcessCombo

/* ----- ProcessMessage ----- */
Bool_t
XSReactionDlg::ProcessMessage(Long_t msg, Long_t param1, Long_t param2)
{
   printf("Message = %d (%d)\n", GET_MSG(msg), GET_SUBMSG(msg));

   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               ProcessButton(param1,param2);
               break;

            case kCM_COMBOBOX:
               ProcessCombo(param1,param2);
               break;

            default:
               break;
         }
      case kC_TEXTENTRY:
         switch (GET_SUBMSG(msg)) {
            case kTE_ENTER:
               ElementEntryChanged();
               break;
            default:
               break;
         }
      default:
//         printf("Message = %d (%d)\n", GET_MSG(msg), GET_SUBMSG(msg));
         break;
   }
   return kTRUE;
} // ProcessMessage
