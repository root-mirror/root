// Author: Roel Aaij 21/07/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGCanvas.h"
#include "TGWindow.h"
#include "TGResourcePool.h"
#include <iostream>
#include <iomanip>
#include "TGWidget.h"
#include "TRandom3.h"
#include "TVirtualTableInterface.h"
#include "TGTable.h"
#include "TGTableCell.h"
#include "TGTableHeader.h"
#include "TObjArray.h"
#include "TGTableContainer.h"
#include "TGScrollBar.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGLabel.h"
#include "TColor.h"

ClassImp(TGTable);
ClassImp(TTableRange);

////////////////////////////////////////////////////////////////////////////////
/// Create an array to hold a bunch of numbers

/* Begin_Html
<center><h2>TGTable</h2></center>
<br><br>
TGTable implements a table widget to display data in rows and
columns. The data is supplied by a TVirtualTableInterface.
<br><br>
The table is a TGCanvas to make use of already available viewport
functionality and drawing optimizations.
<br><br>
The top left cell in a table has coordinates (0,0)
<br><br>
A TObjArray is used internally to ensure little overhead and fast
acces to cells.
<br><br>
If the data source has more rows than the default 50 rows of cells in
memory, buttons at the bottom of the table can be used to load the
next or previous chunk of data.
<br><br>
At the top of the table, a frame is visible that shows the coordinates
of the top left cell currently in memmory in row,column. The amount of
rows and columns is also shown in rows x columns. These values can be
edited to move to a different area of the data source or to resize the
table. Tab will switch between the enties, return will move to the
currently entered range and resize the table if needed. Clicking the
goto button has the same effect.
<br><br>
A TGTable is created by first creating an appropriate
TVirtualTableInterface from the data that needs visualization and
then creating the TGTable using this interface.
<br><br>
A simple macro to use a TGTable with a TGSimpleTableInterface:
End_Html
~~~
{
   Int_t i = 0, j = 0;
   UInt_t nrows = 6, ncolumns = 5;
   Double_t** data = new Double_t*[nrows];
   for (i = 0; i < nrows; i++) {
      data[i] = new Double_t[ncolumns];
      for (j = 0; j < ncolumns; j++) {
         data[i][j] = 10 * i + j;
      }
   }

   // Create a main frame to contain the table
   TGMainFrame* mainframe = new TGMainFrame(0, 400, 200);
   mainframe->SetCleanup(kDeepCleanup) ;

   // Create an interface
   TGSimpleTableInterface *iface = new TGSimpleTableInterface(data, 6, 5);

   // Create the table
   TGTable *table = new TGTable(mainframe, 999, iface);

   // Add the table to the main frame
   mainframe->AddFrame(table, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   //Update data
   data[5][1] = 3.01;
   //update the table view
   table->Update();

   // Layout and map the main frame
   mainframe->SetWindowName("Tree Table Test") ;
   mainframe->MapSubwindows() ;
   mainframe->Layout();
   mainframe->Resize() ;
   mainframe->MapWindow() ;

   return mainframe;
}
~~~
Begin_Html

It is also possible to visualise data from a tree. A simple macro
showing the use of a TTreeTableInterface follows.
End_Html
~~~
{
   // Open a root file.
   TFile *file = new TFile("$ROOTSYS/tutorials/hsimple.root");
   // Load a tree from the file
   TNtuple *ntuple = (TNtuple *)file->Get("ntuple");

   // Create an interface
   TTreeTableInterface *iface = new TTreeTableInterface(ntuple);

   // Create a main frame to contain the table
   TGMainFrame* mainframe = new TGMainFrame(0, 400, 200);
   mainframe->SetCleanup(kDeepCleanup) ;

   // Create the table
   TGTable *table = new TGTable(mainframe, 999, iface, 10, 6);

   // Add the table to the main frame
   mainframe->AddFrame(table, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   // Set a selection
   iface->SetSelection("px > 0.");
   // Add a column
   iface->AddColumn("(px+py)/(px-py)", 0);
   //update the table view
   table->Update();

   // Layout and map the main frame
   mainframe->SetWindowName("Tree Table Test") ;
   mainframe->MapSubwindows() ;
   mainframe->Layout();
   mainframe->Resize() ;
   mainframe->MapWindow() ;

   return mainframe;
}
~~~
*/

// const TGGC *TGTable::fgDefaultSelectGC = 0;
// const TGGC *TGTable::fgDefaultBckgndGC = 0;
// const Int_t TGTable::fgDefaultTMode = kTextLeft | kTextTop;

// TList *TGTable::fgEditList = 0 ;

////////////////////////////////////////////////////////////////////////////////
/// TGTable constuctor.

TGTable::TGTable(const TGWindow *p, Int_t id, TVirtualTableInterface *interface,
                 UInt_t nrows, UInt_t ncolumns)
   : TGCompositeFrame(p, 500, 500, kVerticalFrame), TGWidget(id), fRows(0),
     fRowHeaders(0), fColumnHeaders(0), fReadOnly(kFALSE), fSelectColor(0),
     fTMode(0), fAllData(kFALSE), fTableFrame(0), fCanvas(0), fCellWidth(80),
     fCellHeight(25), fInterface(interface)
{
   fCurrentRange = new TTableRange();
   fDataRange = new TTableRange();
   fGotoRange = new TTableRange();
   TGLayoutHints *hints = 0;
   fCellHintsList = new TList(hints);
   fRHdrHintsList = new TList(hints);
   fCHdrHintsList = new TList(hints);
   fMainHintsList = new TList(hints);

   // To be done: GetBackground colors for .rootrc
   SetBackgroundColor(GetWhitePixel());
   fEvenRowBackground = TColor::RGB2Pixel(204, 255, 204);
   fOddRowBackground  = TColor::RGB2Pixel(255, 255, 255);
   fHeaderBackground  = TColor::RGB2Pixel(204, 204, 255);

   fCurrentRange->fXbr = ncolumns;
   fCurrentRange->fYbr = nrows;

   Init();

   if(fInterface) SetInterface(fInterface, nrows, ncolumns);
   SetWindowName();
//    MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// TGTable destructor.

TGTable::~TGTable()
{
   // delete all cells in a good way
   UInt_t i = 0, j = 0;
   for (i = 0; i < GetNTableRows(); i++) {
      for (j = 0; j < GetNTableColumns(); j++) {
         delete GetCell(i,j);
      }
      delete fRows->At(i);
   }
   delete fRows;
   delete fRowHeaders;
   delete fColumnHeaders;

   delete fCurrentRange;
   delete fDataRange;
   delete fGotoRange;

   fCellHintsList->Delete();
   delete fCellHintsList;
   delete fRHdrHintsList;
   delete fCHdrHintsList;

   fMainHintsList->Delete();
   delete fMainHintsList;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialise the TGTable.

void TGTable::Init()
{
   UInt_t nrows = GetNTableRows();
   UInt_t ncolumns = GetNTableColumns();

   // Main layout frames
   fTopFrame = new TGHorizontalFrame(this, fWidth, fCellHeight);
   fTopExtraFrame = new TGHorizontalFrame(fTopFrame, fWidth - fCellWidth,
                                          fCellHeight);
   TGString *str = new TGString();
   *str += GetNTableRows();
   *str += "x";
   *str += GetNTableColumns();
   *str += " Table";
   fTableHeader = new TGTableHeader(fTopFrame, this, str, 0,
                                    kTableHeader);

   fBottomFrame = new TGHorizontalFrame(this, fWidth, fHeight - fCellHeight);
   fRHdrFrame = new TGTableHeaderFrame(fBottomFrame, this, fCellWidth,
                                       fHeight - fCellHeight, kRowHeader);
   fCHdrFrame = new TGTableHeaderFrame(fTopExtraFrame, this, fWidth - fCellWidth,
                                       fCellHeight, kColumnHeader);

   // Frame for the buttons at the bottom
   fButtonFrame = new TGHorizontalFrame(this, 200, 50);
   fNextButton = new TGTextButton(fButtonFrame, "Next", WidgetId() + 2000);
   fPrevButton = new TGTextButton(fButtonFrame, "Previous", WidgetId() + 2001);
   fUpdateButton = new TGTextButton(fButtonFrame, "Update", WidgetId() + 2002);

   fCanvas = new TGCanvas(fBottomFrame, ncolumns * fCellWidth,
                          nrows * fCellHeight, 0);
   fTableFrame = new TGTableFrame(fCanvas->GetViewPort(), nrows, ncolumns);
   fTableFrame->SetCanvas(fCanvas);
   fCanvas->SetContainer(fTableFrame->GetFrame());

   // Frame to display range info and goto button.
   fRangeFrame = new TGHorizontalFrame(this, 450, 50);
   fFirstCellLabel = new TGLabel(fRangeFrame, "Top left cell in range:");
   fRangeLabel = new TGLabel(fRangeFrame, "Range:");
   fFirstCellEntry = new TGTextEntry(fRangeFrame, "0,0", WidgetId() + 2050);
   fFirstCellEntry->SetWidth(100);
   fFirstCellEntry->SetAlignment(kTextRight);
   fFirstCellEntry->Connect("TextChanged(const char *)", "TGTable", this,
                            "UserRangeChange()");
   fFirstCellEntry->Connect("ReturnPressed()", "TGTable", this, "Goto()");

   TString range;
   range += GetNTableRows();
   range += "x";
   range += GetNTableColumns();
   fRangeEntry = new TGTextEntry(range, fRangeFrame, WidgetId() + 2051);
   fRangeEntry->SetWidth(100);
   fRangeEntry->SetAlignment(kTextRight);
   fRangeEntry->Connect("TextChanged(const char *)", "TGTable", this,
                        "UserRangeChange()");
   fRangeEntry->Connect("ReturnPressed()", "TGTable", this, "Goto()");
   fRangeEntry->Connect("TabPressed()", "TGTextEntry", fFirstCellEntry,
                            "SetFocus()");
   fFirstCellEntry->Connect("TabPressed()", "TGTextEntry", fRangeEntry,
                            "SetFocus()");

   fGotoRange->fXbr = GetNTableRows();
   fGotoRange->fYbr = GetNTableColumns();
   fGotoButton = new TGTextButton(fRangeFrame, "Goto", WidgetId() + 2003);
   fGotoButton->SetState(kButtonDisabled);

   // Set frame backgrounds
   fCHdrFrame->SetBackgroundColor(fBackground);
   fRHdrFrame->SetBackgroundColor(fBackground);
   fRangeFrame->SetBackgroundColor(fBackground);
   fTopFrame->SetBackgroundColor(fBackground);
   fTopExtraFrame->SetBackgroundColor(fBackground);
   fBottomFrame->SetBackgroundColor(fBackground);
   fButtonFrame->SetBackgroundColor(fBackground);
   fFirstCellLabel->SetBackgroundColor(fBackground);
   fRangeLabel->SetBackgroundColor(fBackground);

   // Create the cells needed
   UInt_t i = 0, j = 0;
   TGString *label = 0;
   fRowHeaders = new TObjArray(nrows);
   for(i = 0; i < nrows; i++) {
      TGTableHeader *hdr = new TGTableHeader(fRHdrFrame, this,
                                             label, i, kRowHeader);
      fRowHeaders->AddAt(hdr, i);
   }
   fColumnHeaders = new TObjArray(ncolumns);
   for(i = 0; i < ncolumns; i++) {
      TGTableHeader *hdr = new TGTableHeader(fCHdrFrame, this,
                                             label, i, kColumnHeader);
      fColumnHeaders->AddAt(hdr, i);
   }

   TGTableCell *cell = 0;
   TObjArray *row = 0;
   fRows = new TObjArray(nrows);
   for (i = 0; i < nrows; i++) {
      row = new TObjArray(ncolumns);
      fRows->AddAt(row, i);
      for (j = 0; j < ncolumns; j++) {
         cell = new TGTableCell(fCanvas->GetContainer(), this, label, i, j);
         row->AddAt(cell, j);
      }
   }

   // Check if the table covers all the data
   if ((GetNDataColumns() >= GetNTableColumns()) &&
       (GetNDataRows() >= GetNTableRows())) {
      fAllData = kTRUE;
   } else {
      fAllData = kFALSE;
   }

   TGLayoutHints *lhints = 0;

   // Add cells and headers to layout frames
   for (i = 0; i < nrows; i++) {
      lhints = new TGLayoutHints(kLHintsLeft | kLHintsTop);
      fRHdrHintsList->Add(lhints);
      fRHdrFrame->AddFrame(GetRowHeader(i), lhints);
      for (j = 0; j < ncolumns; j++) {
         if (i == 0) {
            lhints = new TGLayoutHints(kLHintsLeft | kLHintsTop);
            fCHdrHintsList->Add(lhints);
            fCHdrFrame->AddFrame(GetColumnHeader(j), lhints);
         }
         lhints = new TGLayoutHints(kLHintsLeft | kLHintsTop);
         fCellHintsList->Add(lhints);
         fCanvas->AddFrame(GetCell(i,j), lhints);
      }
   }

   // Add frames to the range frame
   lhints = new TGLayoutHints(kLHintsRight | kLHintsCenterY, 3, 30, 4, 4);
   fRangeFrame->AddFrame(fGotoButton, lhints);
   lhints = new TGLayoutHints(kLHintsRight |kLHintsCenterY, 3, 3, 4, 4);
   fRangeFrame->AddFrame(fRangeEntry, lhints);
   lhints = new TGLayoutHints(kLHintsRight |kLHintsCenterY, 3, 3, 4, 4);
   fRangeFrame->AddFrame(fRangeLabel, lhints);
   lhints = new TGLayoutHints(kLHintsRight |kLHintsCenterY, 3, 3, 4, 4);
   fRangeFrame->AddFrame(fFirstCellEntry, lhints);
   lhints = new TGLayoutHints(kLHintsRight |kLHintsCenterY, 3, 3, 4, 4);
   fRangeFrame->AddFrame(fFirstCellLabel, lhints);
   lhints = new TGLayoutHints(kLHintsRight |kLHintsTop);
   fRangeFrame->Resize();
   // Range frame size = 448
   AddFrame(fRangeFrame, lhints);

   // Add table to the main composite frame
   lhints = new TGLayoutHints(kLHintsLeft |kLHintsTop);
   fTopFrame->AddFrame(fTableHeader, lhints);
   lhints = new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsTop);
   fTopExtraFrame->AddFrame(fCHdrFrame, lhints);
   lhints = new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsTop);
   fTopFrame->AddFrame(fTopExtraFrame, lhints);
   lhints = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY);
   fBottomFrame->AddFrame(fRHdrFrame, lhints);
   lhints =  new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX
                               | kLHintsExpandY);
   fBottomFrame->AddFrame(fCanvas, lhints);

   // Add buttons to button frame
   lhints = new TGLayoutHints(kLHintsRight | kLHintsCenterY, 3, 30, 4, 4);
   fButtonFrame->AddFrame(fNextButton, lhints);
   lhints = new TGLayoutHints(kLHintsRight | kLHintsCenterY, 3, 3, 4, 4);
   fButtonFrame->AddFrame(fPrevButton, lhints);
   lhints = new TGLayoutHints(kLHintsRight | kLHintsCenterY, 3, 30, 4, 4);
   fButtonFrame->AddFrame(fUpdateButton, lhints);
   fButtonFrame->Resize();
   fButtonFrame->ChangeOptions(fButtonFrame->GetOptions() | kFixedWidth);

   lhints = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX);
   AddFrame(fTopFrame, lhints);
   lhints = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX |
                              kLHintsExpandY);
   AddFrame(fBottomFrame, lhints);
   lhints = new TGLayoutHints(kLHintsExpandX | kLHintsTop);
   AddFrame(fButtonFrame, lhints);

   // Setup scrolling for the headers
   TGScrollBar *sbar= fCanvas->GetVScrollbar();
   sbar->Connect("PositionChanged(Int_t)", "TGTable", this, "ScrollRHeaders(Int_t)");
   sbar = fCanvas->GetHScrollbar();
   sbar->Connect("PositionChanged(Int_t)", "TGTable", this, "ScrollCHeaders(Int_t)");

   // Connections for buttons
   fUpdateButton->Connect("Clicked()", "TGTable", this, "Update()");
   fNextButton->Connect("Clicked()", "TGTable", this, "NextChunk()");
   fPrevButton->Connect("Clicked()", "TGTable", this, "PreviousChunk()");
   fGotoButton->Connect("Clicked()", "TGTable", this, "Goto()");

//    MapSubwindows();
//    Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the TGTable.

void TGTable::DoRedraw()
{
   MapSubwindows();
   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Expand a TGTable by nrows and ncolumns.

void TGTable::Expand(UInt_t nrows, UInt_t ncolumns)
{
   ExpandRows(nrows);
   ExpandColumns(ncolumns);
}

////////////////////////////////////////////////////////////////////////////////
/// Expand the columns of a TGTable by ncolumns.

void TGTable::ExpandColumns(UInt_t ncolumns)
{
   UInt_t i = 0, j = 0;
   TGString *label = 0;

   UInt_t ntrows = GetNTableRows();
   UInt_t ntcolumns = GetNTableColumns();

   fColumnHeaders->Expand(ntcolumns + ncolumns);

   for (i = 0; i < ncolumns; i++) {
      TGTableHeader *header = new TGTableHeader(fCHdrFrame, this, label,
                                                ntcolumns + i,
                                                kColumnHeader);
      fColumnHeaders->AddAt(header, ntcolumns + i);
   }

   for (i = 0; i < ntrows; i++) {
      GetRow(i)->Expand(ntcolumns + ncolumns);
      for (j = 0; j < ncolumns; j++) {
         TGTableCell *cell = new TGTableCell(fCanvas->GetContainer(), this, label, i,
                                             ntcolumns + j);
         if (GetRow(i)) GetRow(i)->AddAt(cell, ntcolumns + j);
      }
   }

   fCurrentRange->fXbr += ncolumns;

   if ((GetNDataColumns() == GetNTableColumns()) &&
       (GetNDataRows() == GetNTableRows())) {
      fAllData = kTRUE;
   } else {
      fAllData = kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Expand the rows of a TGTable by nrows.

void TGTable::ExpandRows(UInt_t nrows)
{
   UInt_t i = 0, j = 0;

   UInt_t ntrows = GetNTableRows();
   UInt_t ntcolumns = GetNTableColumns();

   fRows->Expand(ntrows + nrows);
   fRowHeaders->Expand(ntrows + nrows);
   for (i = 0; i < nrows; i++) {
      TObjArray *row = new TObjArray(ntcolumns);
      fRows->AddAt(row, ntrows + i);
      TGString *label = 0;
      TGTableHeader *header = new TGTableHeader(fRHdrFrame, this, label,
                                                ntrows + i, kRowHeader);
      fRowHeaders->AddAt(header, ntrows + i);
      for (j = 0; j < ntcolumns ; j++) {
         TGTableCell *cell = new TGTableCell(fCanvas->GetContainer(), this, label,
                                             ntrows + i, j);
         if (GetRow(ntrows + i)) GetRow(ntrows + i)->AddAt(cell, j);
      }
   }

   fCurrentRange->fYbr += nrows;

   if ((GetNDataColumns() == GetNTableColumns()) &&
       (GetNDataRows() == GetNTableRows())) {
      fAllData = kTRUE;
   } else {
      fAllData = kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the current width of the column header frame.

UInt_t TGTable::GetCHdrWidth() const
{
   Int_t ncolumns = GetNTableColumns();
   UInt_t width = 0;
   for (Int_t i = 0; i < ncolumns; i++) {
      if (GetColumnHeader(i)) width += GetColumnHeader(i)->GetWidth();
   }
   return width;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the current height of the row header frame.

UInt_t TGTable::GetRHdrHeight() const
{
   Int_t nrows = GetNTableRows();
   UInt_t height = 0;
   for (Int_t i = 0; i < nrows; i++) {
      if (GetRowHeader(i)) height += GetRowHeader(i)->GetHeight();
   }
   return height;
}

////////////////////////////////////////////////////////////////////////////////
/// Shrink the TGTable by nrows and ncolumns.

void TGTable::Shrink(UInt_t nrows, UInt_t ncolumns)
{
   ShrinkRows(nrows);
   ShrinkColumns(ncolumns);
}

////////////////////////////////////////////////////////////////////////////////
/// Shrink the columns of the TGTable by ncolumns.

void TGTable::ShrinkColumns(UInt_t ncolumns)
{
   UInt_t i = 0, j = 0, k = 0;

   if(GetNTableColumns() - ncolumns < 1) {
      Info("TGTable::ShrinkColumns", "Cannot shrink smaller than 1"
                                     " column, adjusting");
      ncolumns = GetNTableColumns() - 1;
   }

   UInt_t ntrows = GetNTableRows();
   UInt_t ntcolumns = GetNTableColumns();

   TGTableCell *cell = 0;

   //Destroy windows

   for (i = 0; i < ntrows; i++) {
      for (j = 0; j < ncolumns; j++) {
         k = ntcolumns - ncolumns + j;
         if (GetRow(i)) {
            cell = (TGTableCell *)GetRow(i)->RemoveAt(k);
            if (cell) {
               cell->DestroyWindow();
               delete cell;
            }
         }
      }
      GetRow(i)->Expand(ntcolumns - ncolumns);
   }

   TGTableHeader *hdr = 0;
   for (j = 0; j < ncolumns; j++) {
      hdr = (TGTableHeader *)fColumnHeaders->RemoveAt(ntcolumns - ncolumns + j);
      hdr->DestroyWindow();
      delete hdr;
   }
   fColumnHeaders->Expand(ntcolumns - ncolumns);

   fCurrentRange->fXbr -= ncolumns;


   if ((GetNDataColumns() == GetNTableColumns()) &&
       (GetNDataRows() == GetNTableRows())) {
      fAllData = kTRUE;
   } else {
      fAllData = kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Shrink the rows of the TGTable by nrows.

void TGTable::ShrinkRows(UInt_t nrows)
{
   UInt_t i = 0 , j = 0;

   if(GetNTableRows() - nrows < 1) {
      Info("TGTable::ShrinkRows", "Cannot shrink smaller than 1 row, adjusting");
      nrows = GetNTableRows() - 1;
   }

   UInt_t ntrows = GetNTableRows();
   UInt_t ntcolumns = GetNTableColumns();

   TObjArray *row = 0;
   TGTableCell *cell = 0;
   TGTableHeader *hdr = 0;

   for (i = 0; i < nrows; i++) {
      for (j = 0; j < ntcolumns ; j++) {
         if (GetRow(ntrows - nrows + i)) {
            cell = (TGTableCell *)GetRow(ntrows - nrows + i)->RemoveAt(j);
            if (cell) {
               cell->DestroyWindow();
               delete cell;
            }
         }
      }
      row = (TObjArray *)fRows->RemoveAt(ntrows - nrows + i);
      delete row;
      hdr = (TGTableHeader *)fRowHeaders->RemoveAt(ntrows - nrows + i);
      hdr->DestroyWindow();
      delete hdr;
   }
   fRows->Expand(ntrows - nrows);
   fRowHeaders->Expand(ntrows - nrows);

   fCurrentRange->fYbr -= nrows;

   if ((GetNDataColumns() == GetNTableColumns()) &&
       (GetNDataRows() == GetNTableRows())) {
      fAllData = kTRUE;
   } else {
      fAllData = kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update the labels of the headers of the given type

void TGTable::UpdateHeaders(EHeaderType type)
{
   UInt_t max = 0, i = 0, d = 0;
   if(type == kColumnHeader) {
      max = GetNTableColumns();
      for (i = 0; i < max; i++) {
         d = fCurrentRange->fXtl + i;
         if (GetColumnHeader(i) && fInterface->GetColumnHeader(d))
            GetColumnHeader(i)->SetLabel(fInterface->GetColumnHeader(d));
      }
   } else if (type == kRowHeader) {
      max = GetNTableRows();
      for (i = 0; i < max; i++) {
         d = fCurrentRange->fYtl + i;
         if (GetRowHeader(i) && fInterface->GetRowHeader(d))
            GetRowHeader(i)->SetLabel(fInterface->GetRowHeader(d));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the interface that the TGTable uses to interface.

void TGTable::SetInterface(TVirtualTableInterface *interface,
                           UInt_t nrows, UInt_t ncolumns)
{
   fInterface = interface;

   // Set up ranges

   fDataRange->fXtl = 0;
   fDataRange->fYtl = 0;
   fDataRange->fXbr = fInterface->GetNColumns();
   fDataRange->fYbr = fInterface->GetNRows();

   UInt_t x = 0, y = 0;
   if (fDataRange->fXbr < ncolumns) {
      x = fDataRange->fXbr;
   } else {
      x = ncolumns;
   }

   if (fDataRange->fYbr < nrows) {
      y = fDataRange->fYbr;
   } else {
      y = nrows;
   }

   GotoTableRange(0, 0, x, y);

   if ((GetNDataColumns() == GetNTableColumns()) &&
       (GetNDataRows() == GetNTableRows())) {
      fAllData = kTRUE;
   } else {
      fAllData = kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Resize the table to newnrows and newncolumns and add all the frames to
/// their parent frames.

void TGTable::ResizeTable(UInt_t newnrows, UInt_t newncolumns)
{
   UInt_t oldnrows = GetNTableRows();
   UInt_t oldncolumns = GetNTableColumns();

   Int_t i = 0, j = 0;

   TGCompositeFrame *container = (TGCompositeFrame *)fCanvas->GetContainer();

   if (newnrows != oldnrows){
      if (newnrows > oldnrows) {
         ExpandRows(newnrows - oldnrows);
      } else {
         ShrinkRows(oldnrows - newnrows);
      }
   }

   if (newncolumns != oldncolumns){
      if (newncolumns > oldncolumns) {
         ExpandColumns(newncolumns - oldncolumns);
      } else {
         ShrinkColumns(oldncolumns - newncolumns);
      }
   }

   // Update the layoutmanager and add the frames.
   if ((newncolumns != oldncolumns) || (newnrows != oldnrows)) {
      container->RemoveAll();
      fCellHintsList->Delete();

      fRHdrFrame->RemoveAll();
      fRHdrHintsList->Delete();

      fCHdrFrame->RemoveAll();
      fCHdrHintsList->Delete();

      container->SetLayoutManager(new TGMatrixLayout(container,
                                                     newnrows, newncolumns));
      // Add frames to layout frames
      TGLayoutHints *lhints = 0;
      for (i = 0; i < (Int_t)newnrows; i++) {
         lhints = new TGLayoutHints(kLHintsLeft | kLHintsTop);
         fRHdrHintsList->Add(lhints);
         fRHdrFrame->AddFrame(GetRowHeader(i), lhints);
         for (j = 0; j < (Int_t)newncolumns; j++) {
            if (i == 0) {
               lhints = new TGLayoutHints(kLHintsLeft | kLHintsTop);
               fCHdrHintsList->Add(lhints);
               fCHdrFrame->AddFrame(GetColumnHeader(j), lhints);
            }
            lhints = new TGLayoutHints(kLHintsLeft | kLHintsTop);
            fCellHintsList->Add(lhints);
            fCanvas->AddFrame(GetCell(i,j), lhints);
         }
      }
   }
   fCanvas->MapSubwindows();
   fCanvas->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Update the range shown in the range frame.

void TGTable::UpdateRangeFrame()
{
   TString tl, range;

   tl += fCurrentRange->fYtl;
   tl += ",";
   tl += fCurrentRange->fXtl;
   fFirstCellEntry->SetText(tl.Data());

   range += GetNTableRows();
   range += "x";
   range += GetNTableColumns();
   fRangeEntry->SetText(range.Data());

   fGotoButton->SetState(kButtonDisabled);
}

////////////////////////////////////////////////////////////////////////////////
/// Get row. NOTE: Do not delete the TObjArray returned or the cells
/// it contains, they are owned by the TGTable.

TObjArray *TGTable::GetRow(UInt_t row)
{
   return (TObjArray *)fRows->At(row);
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to a TObjArray that contains pointers to all
/// the cells in column. NOTE: The user will have to delete the
/// TObjArray, but do NOT delete the cells it contains, they are
/// owned by the TGTable and will be deleted from the TGTable with
/// undefined consequenses.

TObjArray *TGTable::GetColumn(UInt_t column)
{
   UInt_t nrows = GetNTableRows();

   TObjArray *col = new TObjArray(nrows);
   for(UInt_t ui = 0; ui < nrows; ui++) {
      col->AddAt(GetCell(ui, column), ui);
   }
   return col;
}

// //______________________________________________________________________________
// void TGTable::Select(TGTableCell *celltl, TGTableCell *cellbr)
// {
// }

// //______________________________________________________________________________
// void TGTable::Select(UInt_t xcelltl, UInt_t ycelltl, UInt_t xcell2, UInt_t ycell2)
// {
// }

// //______________________________________________________________________________
// void TGTable::SelectAll()
// {
// }

// //______________________________________________________________________________
// void TGTable::SelectRow(TGTableCell *cell)
// {
// }

// //______________________________________________________________________________
// void TGTable::SelectRow(UInt_t row)
// {
// }

// //______________________________________________________________________________
// void TGTable::SelectRows(UInt_t row, UInt_t nrows)
// {
// }

// //______________________________________________________________________________
// void TGTable::SelectColumn(TGTableCell *cell)
// {
// }

// //______________________________________________________________________________
// void TGTable::SelectColumn(UInt_t column)
// {
// }

// //______________________________________________________________________________
// void TGTable::SelectColumns(UInt_t column, UInt_t ncolumns)
// {
// }

// //______________________________________________________________________________
// void TGTable::SetBckgndGC(TGGC *gc)
// {
// }

// //______________________________________________________________________________
// void TGTable::SetSelectGC(TGGC *gc)
// {
// }

// //______________________________________________________________________________
// void TGTable::SetTextJustify(Int_t tmode)
// {
// }

////////////////////////////////////////////////////////////////////////////////
/// Const version of GetCell().

const TGTableCell* TGTable::GetCell(UInt_t i, UInt_t j) const
{
   return const_cast<TGTable *>(this)->GetCell(i, j);
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the TGTableCell at position i,j.

TGTableCell* TGTable::GetCell(UInt_t i, UInt_t j)
{
   TObjArray *row = (TObjArray *)fRows->At(i);
   if(row) {
      TGTableCell *cell = (TGTableCell *)row->At(j);
      return cell;
   } else {
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Const version of FindCell().

const TGTableCell* TGTable::FindCell(TGString label) const
{
   return const_cast<TGTable *>(this)->FindCell(label);
}

////////////////////////////////////////////////////////////////////////////////
/// Find the TGTableCell with label.

TGTableCell* TGTable::FindCell(TGString label)
{
   TObjArray *row = 0;
   TGTableCell *cell = 0;
   UInt_t i = 0, j = 0;
   //continue here
   UInt_t nrows = GetNTableRows();
   UInt_t ncolumns = GetNTableColumns();
   for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncolumns; j++) {
         row = (TObjArray *)fRows->At(j);
         cell = (TGTableCell *)row->At(i);
         if (*(cell->GetLabel()) == label) {
            return cell;
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Show the contents of the TGTable in stdout.

void TGTable::Show()
{
   TGTableCell *cell = 0;
   TGTableHeader *hdr = 0;
   UInt_t i = 0, j = 0;
   UInt_t nrows = GetNTableRows();
   UInt_t ncolumns = GetNTableColumns();

   // save actual formatting flags
   std::ios_base::fmtflags org_flags = std::cout.flags();

   for (j = 0; j < ncolumns + 1; j++) {
      if (j == 0) {
         hdr = fTableHeader;
         if (hdr) std::cout << " " << std::setw(12) << std::right
                            << hdr->GetLabel()->GetString() << " ";
      } else {
         hdr = GetColumnHeader(j - 1);
         if (hdr) std::cout << " " << std::setw(12) << std::right
                            << hdr->GetLabel()->GetString() << " ";
      }
   }
   std::cout << std::endl;

   for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncolumns + 1; j++) {
         if (j == 0) {
            hdr = GetRowHeader(i);
            if (hdr) std::cout << " " << std::setw(12) << std::right
                               << hdr->GetLabel()->GetString() << " ";
         } else {
            cell = GetCell(i, j - 1);
            if (cell) std::cout << " " << std::setw(12) << std::right
                                << cell->GetLabel()->GetString() << " ";
         }
      }
      std::cout << std::endl;
   }
   // restore original formatting flags
   std::cout.flags(org_flags);
}

// //______________________________________________________________________________
// void TGTable::InsertRowBefore(UInt_t row, UInt_t nrows)
// {
// }

// //______________________________________________________________________________
// void TGTable::InsertRowBefore(TGString label, UInt_t nrows)
// {
// }

// //______________________________________________________________________________
// void TGTable::InsertRowAfter(UInt_t row, UInt_t nrows)
// {
// }

// //______________________________________________________________________________
// void TGTable::InsertRowAfter(TGString label, UInt_t nrows)
// {
// }

// //______________________________________________________________________________
// void TGTable::InsertRowAt(UInt_t row, UInt_t nrows)
// {
// }

// //______________________________________________________________________________
// void TGTable::InsertRowAt(TGString label, UInt_t nrows)
// {
// }

// //______________________________________________________________________________
// void TGTable::InsertColumnBefore(UInt_t column, UInt_t ncolumns)
// {
// }

// //______________________________________________________________________________
// void TGTable::InsertColumnBefore(TGString label, UInt_t ncolumns)
// {
// }

// //______________________________________________________________________________
// void TGTable::InsertColumnAfter(UInt_t column, UInt_t ncolumns)
// {
// }

// //______________________________________________________________________________
// void TGTable::InsertColumnAfter(TGString label, UInt_t ncolumns)
// {
// }

// //______________________________________________________________________________
// void TGTable::InsertColumnAt(UInt_t column, UInt_t ncolumns)
// {
// }

// //______________________________________________________________________________
// void TGTable::InsertColumnAt(TGString label, UInt_t ncolumns)
// {
// }

// //______________________________________________________________________________
// void TGTable::RemoveRows(UInt_t row, UInt_t nrows)
// {
// }

// //______________________________________________________________________________
// void TGTable::RemoveColumns(UInt_t column, UInt_t ncolumns)
// {
// }

////////////////////////////////////////////////////////////////////////////////
/// Update and layout the visible part of the TGTable.

void TGTable::UpdateView()
{
   UInt_t nrows = GetNTableRows();
   UInt_t ncolumns = GetNTableColumns();

   TGString *str = new TGString();
   *str += nrows;
   *str += "x";
   *str += ncolumns;
   *str += " Table";
   fTableHeader->SetLabel(str->GetString());
   delete str;

   UpdateHeaders(kRowHeader);
   UpdateHeaders(kColumnHeader);

   UInt_t i = 0, j = 0;
   UInt_t k = 0, l = 0;

   TGTableCell * cell = 0;
   for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncolumns; j++) {
         cell = GetCell(i,j);
         k = fCurrentRange->fYtl + i;
         l = fCurrentRange->fXtl + j;

         const char *label = fInterface->GetValueAsString(k,l);
         if(cell) cell->SetLabel(label);
      }
   }

   MapSubwindows();
   Layout();
   gClient->NeedRedraw(fTableHeader);
   TGViewPort *vp = fCanvas->GetViewPort();
   fTableFrame->DrawRegion(0, 0, vp->GetWidth(), vp->GetHeight());
   fCHdrFrame->DrawRegion(0, 0, fCHdrFrame->GetWidth(), fCHdrFrame->GetHeight());
   fRHdrFrame->DrawRegion(0, 0, fRHdrFrame->GetWidth(), fRHdrFrame->GetHeight());

   UpdateRangeFrame();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the amount of rows in the table.

UInt_t TGTable::GetNTableRows() const
{
   return fCurrentRange->fYbr - fCurrentRange->fYtl;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the amount of rows in the data source.

UInt_t TGTable::GetNDataRows() const
{
   return fDataRange->fYbr - fDataRange->fYtl;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the amount of columns in the table.

UInt_t TGTable::GetNTableColumns() const
{
   return fCurrentRange->fXbr - fCurrentRange->fXtl;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the amount of columns in the data source.

UInt_t TGTable::GetNDataColumns() const
{
   return fDataRange->fYbr - fDataRange->fYtl;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the amount of cells in the table.

UInt_t TGTable::GetNTableCells() const
{
   return GetNTableRows() * GetNTableColumns();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the amount of cell in the data source.

UInt_t TGTable::GetNDataCells() const
{
   return GetNDataRows() * GetNDataColumns();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current range of the TGTable.

const TTableRange *TGTable::GetCurrentRange() const
{
   return fCurrentRange;
}

////////////////////////////////////////////////////////////////////////////////
/// Const version of GetRowHeader();

const TGTableHeader *TGTable::GetRowHeader(const UInt_t row) const
{
   return const_cast<TGTable *>(this)->GetRowHeader(row);
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the header of row.

TGTableHeader *TGTable::GetRowHeader(const UInt_t row)
{
   return (TGTableHeader *)fRowHeaders->At(row);
}

////////////////////////////////////////////////////////////////////////////////
/// Const version of GetColumnHeader();

const TGTableHeader *TGTable::GetColumnHeader(const UInt_t column) const
{
   return const_cast<TGTable *>(this)->GetColumnHeader(column);
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the header of column.

TGTableHeader *TGTable::GetColumnHeader(const UInt_t column)
{
   return (TGTableHeader *)fColumnHeaders->At(column);
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the table header.

TGTableHeader *TGTable::GetTableHeader()
{
   return fTableHeader;
}

// //______________________________________________________________________________
// const TGGC*  TGTable::GetSelectGC() const
// {
// }

// //______________________________________________________________________________
// const TGGC*  TGTable::GetCellBckgndGC(TGTableCell *cell) const
// {
// }

// //______________________________________________________________________________
// const TGGC*  TGTable::GetCellBckgndGC(UInt_t row, UInt_t column) const
// {
// }

////////////////////////////////////////////////////////////////////////////////
/// Get the background collor for row.

Pixel_t TGTable::GetRowBackground(UInt_t row) const
{
   if (row % 2 == 0) { // Even rows
      return fEvenRowBackground;
   } else {            // Odd rows
      return fOddRowBackground;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the background color of headers.

Pixel_t TGTable::GetHeaderBackground() const
{
   return fHeaderBackground;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the background color for all odd numbered rows.

void TGTable::SetOddRowBackground(Pixel_t pixel)
{
   if(pixel == fOddRowBackground) return;

   fOddRowBackground = pixel;

   UInt_t nrows = GetNTableRows();
   UInt_t ncolumns = GetNTableColumns();
   UInt_t i = 0, j = 0;
   TGTableCell *cell = 0;

   for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncolumns; j++) {
         if (i % 2) {
            cell = GetCell(i,j);
            if (cell) cell->SetBackgroundColor(fOddRowBackground);
         }
      }
   }

   UInt_t width = fCanvas->GetViewPort()->GetWidth();
   UInt_t height = fCanvas->GetViewPort()->GetHeight();
   fTableFrame->DrawRegion(0, 0, width, height);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the background color for all even numbered rows.

void TGTable::SetEvenRowBackground(Pixel_t pixel)
{
   if(pixel == fEvenRowBackground) return;

   fEvenRowBackground = pixel;

   UInt_t nrows = GetNTableRows();
   UInt_t ncolumns = GetNTableColumns();
   UInt_t i = 0, j = 0;
   TGTableCell *cell = 0;

   for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncolumns; j++) {
         if (!(i % 2)) {
            cell = GetCell(i,j);
            if (cell) cell->SetBackgroundColor(fEvenRowBackground);
         }
      }
   }
   UInt_t width = fCanvas->GetViewPort()->GetWidth();
   UInt_t height = fCanvas->GetViewPort()->GetHeight();
   fTableFrame->DrawRegion(0, 0, width, height);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the background color for the headers.

void TGTable::SetHeaderBackground(Pixel_t pixel)
{
   if(pixel == fHeaderBackground) return;

   fHeaderBackground = pixel;

   UInt_t nrows = GetNTableRows();
   UInt_t ncolumns = GetNTableColumns();
   UInt_t i = 0, j = 0;
   TGTableHeader *hdr = 0;

   for (i = 0; i < nrows; i++) {
      hdr = GetRowHeader(i);
      if (hdr) hdr->SetBackgroundColor(fHeaderBackground);
   }
   UInt_t height = fCanvas->GetViewPort()->GetHeight();
   UInt_t width = fTableHeader->GetWidth();
   fRHdrFrame->DrawRegion(0, 0, width, height);

   for (j = 0; j < ncolumns; j++) {
      hdr = GetColumnHeader(j);
      if (hdr) hdr->SetBackgroundColor(fHeaderBackground);
//       gClient->NeedRedraw(hdr);
   }
   width = fCanvas->GetViewPort()->GetWidth();
   height = fTableHeader->GetHeight();
   fCHdrFrame->DrawRegion(0, 0, width, height);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the background color for all rows and headers to their defaults.

void TGTable::SetDefaultColors()
{
   SetEvenRowBackground(TColor::RGB2Pixel(204, 255, 204));
   SetOddRowBackground(TColor::RGB2Pixel(255, 255, 255));
   SetHeaderBackground(TColor::RGB2Pixel(204, 204, 255));
}

////////////////////////////////////////////////////////////////////////////////
/// Move and layout the table to the specified range.

void TGTable::MoveTable(Int_t rows, Int_t columns)
{
   if (fAllData) return;

   Int_t xtl = fCurrentRange->fXtl + columns;
   Int_t ytl = fCurrentRange->fYtl + rows;
   Int_t xbr = fCurrentRange->fXbr + columns;
   Int_t ybr = fCurrentRange->fYbr + rows;

   GotoTableRange(xtl, ytl, xbr, ybr);
}

////////////////////////////////////////////////////////////////////////////////
/// Move and resize the table to the specified range.

void TGTable::GotoTableRange(Int_t xtl,  Int_t ytl, Int_t xbr,  Int_t ybr)
{
   if (fAllData) return;

   if(xtl == xbr || ytl == ybr) {
      Error("TGTable::GotoTableRange","x or y range = 0");
      return;
   }

   Int_t nrows    = TMath::Abs(ybr - ytl);
   Int_t ncolumns = TMath::Abs(xbr - xtl);

   if (xtl > xbr) {
      Info("TGTable::GotoTableRange","Swapping x-range boundries");
      Int_t temp = xtl;
      xtl = xbr;
      xbr = temp;
   }
   if (ytl > ybr) {
      Info("TGTable::GotoTableRange","Swapping y-range boundries");
      Int_t temp = ytl;
      ytl = ybr;
      ybr = temp;
   }

   if((xtl < 0) || (xbr < 0)) {
      Info("TGTable::GotoTableRange", "Column boundry out of bounds, adjusting");
      xtl = 0;
      xbr = ncolumns;
      if (xbr > (Int_t)fDataRange->fXbr) {
         xbr = fDataRange->fXbr;
         ncolumns = TMath::Abs(xbr - xtl);
      }
   }

   if((ytl < 0) || (ybr < 0)) {
      Info("TGTable::GotoTableRange", "Row boundry out of bounds, adjusting");
      ytl = 0;
      ybr = nrows;
      if (ybr > (Int_t)fDataRange->fYbr) {
         ybr = fDataRange->fYbr;
         nrows =  TMath::Abs(ybr - ytl);
      }
   }

   if((xtl > (Int_t)fDataRange->fXbr) || (xbr > (Int_t)fDataRange->fXbr)) {
      Info("TGTable::GotoTableRange", "Left Column boundry out of bounds, "
           "adjusting");
      xbr = fDataRange->fXbr;
      xtl = xbr - ncolumns;
      if (xtl < 0) {
         xtl = 0;
         ncolumns = TMath::Abs(xbr - xtl);
         Info("TGTable::GotoTableRange", "Right column boundry out of"
                                         " bounds, set to 0");
      }
   }
   if ((ytl > (Int_t)fDataRange->fYbr) || (ybr > (Int_t)fDataRange->fYbr)) {
      Info("TGTable::GotoTableRange", "Bottom row boundry out of bounds, "
                                      "adjusting");
      ybr = fDataRange->fYbr;
      ytl = ybr - nrows;
      if (ytl < 0) {
         ytl = 0;
         nrows = ybr - ytl;
         Info("TGTable::GotoTableRange", "Top row boundry out of bounds, "
                                         "set to 0");
      }
   }

   nrows    = TMath::Abs(ybr - ytl);
   ncolumns = TMath::Abs(xbr - xtl);

   // Resize rows and columns if needed
   ResizeTable(nrows, ncolumns);

   fCurrentRange->fXtl = xtl;
   fCurrentRange->fYtl = ytl;
   fCurrentRange->fXbr = xbr;
   fCurrentRange->fYbr = ybr;

   // Update the table view.
   UpdateView();
}

////////////////////////////////////////////////////////////////////////////////
/// Operator for easy cell acces.

TGTableCell *TGTable::operator() (UInt_t row, UInt_t column)
{
   return GetCell(row, column);
}

////////////////////////////////////////////////////////////////////////////////
/// Scroll the column headers horizontally.

void TGTable::ScrollCHeaders(Int_t xpos)
{
   if (!fCHdrFrame) return;

   fCHdrFrame->Move(- xpos, 0);
   fCHdrFrame->Resize();
   fCHdrFrame->DrawRegion(0, 0, fCHdrFrame->GetWidth(),
                          fCHdrFrame->GetHeight());
}

////////////////////////////////////////////////////////////////////////////////
/// Scroll the row headers vertically

void TGTable::ScrollRHeaders(Int_t ypos)
{
   if (!fRHdrFrame) return;

   fRHdrFrame->Move(fRHdrFrame->GetX(), -ypos);
   fRHdrFrame->Resize();
   fRHdrFrame->DrawRegion(0, 0, fRHdrFrame->GetWidth(),
                          fRHdrFrame->GetHeight());
}

////////////////////////////////////////////////////////////////////////////////
/// Move the table to the next chunk of the data set with the same size.

void TGTable::NextChunk()
{
   MoveTable(GetNTableRows(), 0);
   UpdateRangeFrame();
}

////////////////////////////////////////////////////////////////////////////////
/// Move the table to the previous chunk of the data set with the same size.

void TGTable::PreviousChunk()
{
   MoveTable(-1 * (Int_t)GetNTableRows(), 0);
   UpdateRangeFrame();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot used by the Goto button and whenever return is pressed in
/// on of the text entries in the range frame.

void TGTable::Goto()
{
   if (fGotoButton->GetState() == kButtonUp) {
      GotoTableRange(fGotoRange->fXtl, fGotoRange->fYtl,
                     fGotoRange->fXbr, fGotoRange->fYbr);
      UpdateRangeFrame();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot used when the text in one of the range frame text entries changes.

void TGTable::UserRangeChange()
{
   TString topleft(fFirstCellEntry->GetText());
   if(!topleft.Contains(",")) return;

   Int_t pos = topleft.First(',');
   TString itl = topleft(0,pos);
   TString jtl = topleft(pos+1, topleft.Length());

   if (itl.Contains(' ') || itl.Contains('\t') ||
       jtl.Contains(' ') || jtl.Contains('\t')) return;

   if (!itl.IsAlnum() || !jtl.IsAlnum()) return;

   fGotoRange->fXtl = jtl.Atoi();
   fGotoRange->fYtl = itl.Atoi();

   TString range(fRangeEntry->GetText());
   if(!range.Contains("x")) return;

   pos = 0;
   pos = range.First('x');
   TString ir = range(0,pos);
   TString jr = range(pos+1, range.Length());

   if (ir.Contains(' ') || ir.Contains('\t') ||
       jr.Contains(' ') || jr.Contains('\t')) return;
   if (!ir.IsAlnum() || !jr.IsAlnum()) return;

   fGotoRange->fXbr = jtl.Atoi() + jr.Atoi();
   fGotoRange->fYbr = itl.Atoi() + ir.Atoi();

   if (*fGotoRange == *fCurrentRange) {
      fGotoButton->SetState(kButtonDisabled);
   } else {
      fGotoButton->SetState(kButtonUp);
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Update the range of the available data and refresh the current view.

void TGTable::Update()
{
   fDataRange->fXbr = fInterface->GetNColumns();
   fDataRange->fYbr = fInterface->GetNRows();

   GotoTableRange(fCurrentRange->fXtl, fCurrentRange->fYtl,
                  fCurrentRange->fXbr, fCurrentRange->fYbr);

   UpdateView();
}

////////////////////////////////////////////////////////////////////////////////
/// TTableRange constuctor.

TTableRange::TTableRange() : fXtl(0), fYtl(0), fXbr(0), fYbr(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Print the values of a range.

void TTableRange::Print()
{
   std::cout << "Range = (" << fXtl << "," << fYtl << ")->("
             << fXbr << "," << fYbr << ")" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator to determine if 2 ranges are equal

Bool_t TTableRange::operator==(TTableRange &other)
{
   if ((fXtl == other.fXtl) && (fYtl == other.fYtl) &&
       (fXbr == other.fXbr) && (fYbr == other.fYbr)) {
      return kTRUE;
   } else {
      return kFALSE;
   }
}
