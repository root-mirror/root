// @(#)root/treeviewer:$Id$
//Author : Andrei Gheata   16/08/00

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeViewer
#define ROOT_TTreeViewer

////////////////////////////////////////////////////
//                                                //
// TTreeViewer - A GUI oriented tree viewer       //
//                                                //
////////////////////////////////////////////////////

#include "TGFrame.h"

#include "TTree.h"

class TTVLVContainer;
class TTVLVEntry;
class TTVSession;
class TGSelectBox;
class TBranch;
class TContextMenu;
class TList;
class TGPicture;
class TTimer;
class TGLayoutHints;
class TGMenuBar;
class TGPopupMenu;
class TGToolBar;
class TGLabel;
class TGCheckButton;
class TGComboBox;
class TGTextButton;
class TGTextEntry;
class TGDoubleVSlider;
class TGPictureButton;
class TGStatusBar;
class TGCanvas;
class TGListTree;
class TGListTreeItem;
class TGListView;
class TGHProgressBar;
class TGButton;


class TTreeViewer : public TGMainFrame {

friend class TGClient;
friend class TGButton;

public:
   /// Item types used as user data
   enum EListItemType {
      kLTNoType            = 0,
      kLTPackType          = BIT(0),
      kLTTreeType          = BIT(1),
      kLTBranchType        = BIT(2),
      kLTLeafType          = BIT(3),
      kLTActionType        = BIT(4),
      kLTDragType          = BIT(5),
      kLTExpressionType    = BIT(6),
      kLTCutType           = BIT(7)
   };

private:

///@{
   TTree                *fTree;                 ///< Selected tree
   TTVSession           *fSession;              ///< Current tree-viewer session
   const char           *fFilename;             ///< Name of the file containing the tree
   const char           *fSourceFile;           ///< Name of the C++ source file - default treeviewer.C
   TString              fLastOption;            ///< Last graphic option
   TTree                *fMappedTree;           ///< Listed tree
   TBranch              *fMappedBranch;         ///< Listed branch
   Int_t                fDimension;             ///< Histogram dimension
   Bool_t               fVarDraw;               ///< True if an item is double-clicked
   Bool_t               fScanMode;              ///< Flag activated when Scan Box is double-clicked
   TContextMenu         *fContextMenu;          ///< Context menu for tree viewer
   TGSelectBox          *fDialogBox;            ///< Expression editor
   TList                *fTreeList;             ///< List of mapped trees
   Int_t                fTreeIndex;             ///< Index of current tree in list
   const TGPicture      *fPicX;                 ///< Pictures for X expressions
   const TGPicture      *fPicY;                 ///< Pictures for Y expressions
   const TGPicture      *fPicZ;                 ///< Pictures for Z expressions
   const TGPicture      *fPicDraw;              ///< Pictures for Draw buttons
   const TGPicture      *fPicStop;              ///< Pictures for Stop buttons
   const TGPicture      *fPicRefr;              ///< Pictures for Refresh buttons ///<ia
   Cursor_t             fDefaultCursor;         ///< Default cursor
   Cursor_t             fWatchCursor;           ///< Watch cursor
   TTimer               *fTimer;                ///< Tree viewer timer
   Bool_t               fCounting;              ///< True if timer is counting
   Bool_t               fStopMapping;           ///< True if branch don't need remapping
   Bool_t               fEnableCut;             ///< True if cuts are enabled
   Int_t                fNexpressions;          ///< Number of expression widgets
///@}

///@{
/// @name Menu bar, menu bar entries and layouts
   TGLayoutHints        *fMenuBarLayout;
   TGLayoutHints        *fMenuBarItemLayout;
   TGLayoutHints        *fMenuBarHelpLayout;
   TGMenuBar            *fMenuBar;
   TGPopupMenu          *fFileMenu;
   TGPopupMenu          *fEditMenu;
   TGPopupMenu          *fRunMenu;
   TGPopupMenu          *fOptionsMenu;
   TGPopupMenu          *fOptionsGen;
   TGPopupMenu          *fOptions1D;
   TGPopupMenu          *fOptions2D;
   TGPopupMenu          *fHelpMenu;
///@}

///@{
/// @name Toolbar and hints
   TGToolBar            *fToolBar;
   TGLayoutHints        *fBarLayout;
///@}

///@{
/// @name Widgets on the toolbar
   TGLabel              *fBarLbl1;      ///< Label of command text entry
   TGLabel              *fBarLbl2;      ///< Label of option text entry
   TGLabel              *fBarLbl3;      ///< Label of histogram name text entry
   TGCheckButton        *fBarH;         ///< Checked for drawing current histogram with different graphic option
   TGCheckButton        *fBarScan;      ///< Checked for tree scan
   TGCheckButton        *fBarRec;       ///< Command recording toggle
   TGTextEntry          *fBarCommand;   ///< User command entry
   TGTextEntry          *fBarOption;    ///< Histogram drawing option entry
   TGTextEntry          *fBarHist;      ///< Histogram name entry
///@}

///@{
/// @name Frames
   TGHorizontalFrame    *fHf;           ///< Main horizontal frame
   TGDoubleVSlider      *fSlider;       ///< Vertical slider to select processed tree entries;
   TGVerticalFrame      *fV1;           ///< List tree mother
   TGVerticalFrame      *fV2;           ///< List view mother
   TGCompositeFrame     *fTreeHdr;      ///< Header for list tree
   TGCompositeFrame     *fListHdr;      ///< Header for list view
   TGLabel              *fLbl1;         ///< Label for list tree
   TGLabel              *fLbl2;         ///< Label for list view
   TGHorizontalFrame    *fBFrame;       ///< Button frame
   TGHorizontalFrame    *fHpb;          ///< Progress bar frame
   TGHProgressBar       *fProgressBar;  ///< Progress bar
   TGLabel              *fBLbl4;        ///< Label for input list entry
   TGLabel              *fBLbl5;        ///< Label for output list entry
   TGTextEntry          *fBarListIn;    ///< Tree input event list name entry
   TGTextEntry          *fBarListOut;   ///< Pree output event list name entry
   TGPictureButton      *fDRAW;         ///< DRAW button
   TGTextButton         *fSPIDER;       ///< SPIDER button
   TGPictureButton      *fSTOP;         ///< Interrupt current command (not yet)
   TGPictureButton      *fREFR;         ///< REFRESH button  ///<ia
   TGStatusBar          *fStatusBar;    ///< Status bar
   TGComboBox           *fCombo;        ///< Combo box with session records
   TGPictureButton      *fBGFirst;
   TGPictureButton      *fBGPrevious;
   TGPictureButton      *fBGRecord;
   TGPictureButton      *fBGNext;
   TGPictureButton      *fBGLast;
   TGTextButton         *fReset;        ///< clear expression's entries
///@}

///@{
/// @name  ListTree
   TGCanvas             *fTreeView;     ///< ListTree canvas container
   TGListTree           *fLt;           ///< ListTree with file and tree items
///@}

///@{
/// @name ListView
   TGListView           *fListView;     ///< ListView with branches and leaves
   TTVLVContainer       *fLVContainer;  ///< Container for listview

   TList                *fWidgets;      ///< List of widgets to be deleted
///@}


private:

   void          BuildInterface();
   const char   *Cut();
   Int_t         Dimension();
   const char   *EmptyBrackets(const char* name);
   const char   *Ex();
   const char   *Ey();
   const char   *Ez();
   const char   *En(Int_t n);
   void          MapBranch(TBranch *branch, const char *prefix="", TGListTreeItem *parent = 0, Bool_t listIt = kTRUE);
   void          MapOptions(Long_t parm1);
   void          MapTree(TTree *tree, TGListTreeItem *parent = 0, Bool_t listIt = kTRUE);
   void          SetFile();
   const char   *ScanList();
   void          SetParentTree(TGListTreeItem *item);
   void          DoError(int level, const char *location, const char *fmt, va_list va) const;

public:
   TTreeViewer(const char* treeName = 0);
   TTreeViewer(const TTree *tree);
   virtual       ~TTreeViewer();

   void          AppendTree(TTree *tree);
   void          ActivateButtons(Bool_t first, Bool_t previous,
                                 Bool_t next , Bool_t last);
   virtual void  CloseWindow();
   virtual void  Delete(Option_t *) { }                          // *MENU*
   void          DoRefresh();
   void          EditExpression();
   void          Empty();
   void          EmptyAll();                                     // *MENU*
   void          ExecuteCommand(const char* command, Bool_t fast = kFALSE); // *MENU*
   void          ExecuteDraw();
   void          ExecuteSpider();
   TTVLVEntry   *ExpressionItem(Int_t index);
   TList        *ExpressionList();
   const char   *GetGrOpt();
   TTree        *GetTree() {return fTree;}
   Bool_t        HandleTimer(TTimer *timer);
   Bool_t        IsCutEnabled() {return fEnableCut;}
   Bool_t        IsScanRedirected();
   Int_t         MakeSelector(const char* selector = 0);         // *MENU*
   void          Message(const char* msg);
   void          NewExpression();                                // *MENU*
   void          PrintEntries();
   Long64_t      Process(const char* filename, Option_t *option="", Long64_t nentries=TTree::kMaxEntries, Long64_t firstentry=0); // *MENU*
   Bool_t        ProcessMessage(Longptr_t msg, Longptr_t parm1, Longptr_t parm2);
   void          RemoveItem();
   void          RemoveLastRecord();                             // *MENU*
   void          SaveSource(const char* filename="", Option_t *option="");            // *MENU*
   void          SetHistogramTitle(const char *title);
   void          SetCutMode(Bool_t enabled = kTRUE) {fEnableCut = enabled;}
   void          SetCurrentRecord(Long64_t entry);
   void          SetGrOpt(const char *option);
   void          SetNexpressions(Int_t expr);
   void          SetRecordName(const char *name);                // *MENU*
   void          SetScanFileName(const char *name="");           // *MENU*
   void          SetScanMode(Bool_t mode=kTRUE) {fScanMode = mode;}
   void          SetScanRedirect(Bool_t mode);
   void          SetSession(TTVSession *session);
   void          SetUserCode(const char *code, Bool_t autoexec=kTRUE); // *MENU*
   void          SetTree(TTree* tree);
   void          SetTreeName(const char* treeName);              // *MENU*
   Bool_t        SwitchTree(Int_t index);
   void          UpdateCombo();
   void          UpdateRecord(const char *name="new name");      // *MENU*

   ClassDef(TTreeViewer,0)  // A GUI oriented tree viewer
};

#endif
