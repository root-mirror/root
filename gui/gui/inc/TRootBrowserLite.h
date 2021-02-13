// @(#)root/gui:$Id$
// Author: Fons Rademakers   27/02/98

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT_TRootBrowserLite
#define ROOT_TRootBrowserLite

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootBrowserLite                                                     //
//                                                                      //
// This class creates a ROOT object browser (looking like Windows       //
// Explorer). The widgets used are the new native ROOT GUI widgets.     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TBrowserImp.h"
#include "TGFrame.h"

class TGMenuBar;
class TGPopupMenu;
class TGLayoutHints;
class TGStatusBar;
class TGHorizontal3DLine;
class TGToolBar;
class TGButton;
class TGFSComboBox;
class TGLabel;
class TGListView;
class TRootIconBox;
class TGCanvas;
class TGListTree;
class TGListTreeItem;
class TGFileItem;
class TList;
class TGFileContainer;
class TGComboBox;
class TGTextEdit;

class TRootBrowserLite : public TGMainFrame, public TBrowserImp {

friend class TRootIconBox;

private:
   TGMenuBar           *fMenuBar;
   TGToolBar           *fToolBar;
   TGHorizontal3DLine  *fToolBarSep;
   TGVerticalFrame     *fV1;
   TGVerticalFrame     *fV2;
   TGLabel             *fLbl1;
   TGLabel             *fLbl2;
   TGHorizontalFrame   *fHf;
   TGCompositeFrame    *fTreeHdr;
   TGCompositeFrame    *fListHdr;

   TGLayoutHints       *fMenuBarLayout;
   TGLayoutHints       *fMenuBarItemLayout;
   TGLayoutHints       *fMenuBarHelpLayout;
   TGLayoutHints       *fComboLayout;
   TGLayoutHints       *fBarLayout;
   TGComboBox          *fDrawOption;         // drawing option entry
   TGLayoutHints       *fExpandLayout;       //
   Bool_t               fBrowseTextFile;     //
   TString              fTextFileName;

   TList               *fWidgets;
   TList               *fHistory;            // history of browsing
   TObject             *fHistoryCursor;      // current hsitory position
   const TGPicture     *fIconPic;            // icon picture

   void  CreateBrowser(const char *name);
   void  ListTreeHighlight(TGListTreeItem *item);
   void  DeleteListTreeItem(TGListTreeItem *item);
   void  HighlightListLevel();
   void  AddToHistory(TGListTreeItem *item);
   void  IconBoxAction(TObject *obj);
   void  Chdir(TGListTreeItem *item);
   void  DisplayDirectory();
   void  DisplayTotal(Int_t total, Int_t selected);
   void  SetViewMode(Int_t new_mode, Bool_t force = kFALSE);
   void  ToSystemDirectory(const char *dirname);
   void  UpdateDrawOption();
   void  Search();
   void  BrowseTextFile(const char *file);
   void  HideTextEdit();
   void  ShowMacroButtons(Bool_t show = kTRUE);

   Bool_t  HistoryBackward();
   Bool_t  HistoryForward();
   void    ClearHistory();

   TRootBrowserLite(const TRootBrowserLite&) = delete;
   TRootBrowserLite& operator=(const TRootBrowserLite&) = delete;

protected:
   TGPopupMenu         *fFileMenu;
   TGPopupMenu         *fViewMenu;
   TGPopupMenu         *fOptionMenu;
   TGPopupMenu         *fHelpMenu;
   TGPopupMenu         *fSortMenu;
   TGListView          *fListView;
   TRootIconBox        *fIconBox;
   TGCanvas            *fTreeView;
   TGListTree          *fLt;
   TGButton            *fToolBarButton[7];  // same size as gToolBarData[]
   TGFSComboBox        *fFSComboBox;
   TGStatusBar         *fStatusBar;
   TGListTreeItem      *fListLevel;         // current TGListTree level
   Bool_t               fTreeLock;          // true when we want to lock TGListTree
   Int_t                fViewMode;          // current IconBox view mode
   Int_t                fSortMode;          // current IconBox sort mode
   TGTextEdit          *fTextEdit;          // contents of browsed text file

public:
   TRootBrowserLite(TBrowser *b = nullptr, const char *title = "ROOT Browser", UInt_t width = 800, UInt_t height = 500);
   TRootBrowserLite(TBrowser *b, const char *title, Int_t x, Int_t y, UInt_t width, UInt_t height);
   virtual ~TRootBrowserLite();

   void         Add(TObject *obj, const char *name = nullptr, Int_t check = -1) override;
   virtual void AddToBox(TObject *obj, const char *name);
   virtual void AddToTree(TObject *obj, const char *name, Int_t check = -1);

   void         AddCheckBox(TObject *obj, Bool_t check = kFALSE) override;
   void         CheckObjectItem(TObject *obj, Bool_t check = kFALSE) override;
   void         RemoveCheckBox(TObject *obj) override;

   void         BrowseObj(TObject *obj) override;             //*SIGNAL*
   void         ExecuteDefaultAction(TObject *obj) override;  //*SIGNAL*
   virtual void DoubleClicked(TObject *obj);         //*SIGNAL*
   virtual void Checked(TObject *obj, Bool_t check); //*SIGNAL*
   void         CloseTabs() override { }
   void         Iconify() override { }
   void         RecursiveRemove(TObject *obj) override;
   void         Refresh(Bool_t force = kFALSE) override;
   virtual void ResizeBrowser() { }
   virtual void ShowToolBar(Bool_t show = kTRUE);
   virtual void ShowStatusBar(Bool_t show = kTRUE);
   void         Show() override { MapRaised(); }
   virtual void SetDefaults(const char *iconStyle = nullptr, const char *sortBy = nullptr);
   Bool_t       HandleKey(Event_t *event) override;
   void         SetStatusText(const char *txt, Int_t col) override;

   TGListTree      *GetListTree()  const { return fLt; }
   TGFileContainer *GetIconBox()   const;
   TGStatusBar     *GetStatusBar() const { return fStatusBar; }
   TGMenuBar       *GetMenuBar()   const { return fMenuBar; }
   TGToolBar       *GetToolBar()   const { return fToolBar; }
   void             SetDrawOption(Option_t *option = "") override;
   Option_t        *GetDrawOption() const override;
   void             SetSortMode(Int_t new_mode);
   TGMainFrame     *GetMainFrame() const override { return (TGMainFrame *)this; }

   // overridden from TGMainFrame
   void     CloseWindow() override;
   Bool_t   ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) override;
   void     ReallyDelete() override;

   // auxilary (a la private) methods
   void     ExecMacro();
   void     InterruptMacro();

   static TBrowserImp *NewBrowser(TBrowser *b = nullptr, const char *title = "ROOT Browser", UInt_t width = 800, UInt_t height = 500, Option_t *opt="");
   static TBrowserImp *NewBrowser(TBrowser *b, const char *title, Int_t x, Int_t y, UInt_t width, UInt_t height, Option_t *opt="");

   ClassDefOverride(TRootBrowserLite,0)  //ROOT native GUI version of browser
};

#endif
