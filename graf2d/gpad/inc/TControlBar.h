// @(#)root/gpad:$Id$
// Author: Nenad Buncic   20/02/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TControlBar
#define ROOT_TControlBar

#include "TControlBarButton.h"

#ifdef R__LESS_INCLUDES
class TObjArray;
class TControlBarImp;
#else
#include "TControlBarImp.h"
#endif

class TList;

class TControlBar : public TControlBarButton {

friend class  TControlBarImp;

private:
   virtual void   Create();
   void           Initialize(Int_t x, Int_t y);

protected:

   TControlBarImp *fControlBarImp;  ///< system specific implementation
   Int_t           fOrientation;    ///< orientation
   TList          *fButtons;        ///< list of buttons
   Int_t           fNoroc;          ///< number of rows or columns

public:
   enum { kVertical = 1, kHorizontal = 2 };

   TControlBar();
   TControlBar(const char *orientation, const char *title="");
   TControlBar(const char *orientation, const char *title, Int_t x, Int_t y);
   virtual ~TControlBar();

   void            AddButton(TControlBarButton *button);
   void            AddButton(TControlBarButton &button);
   void            AddButton(const char *label, const char *action, const char *hint="", const char *type="button");
   void            AddControlBar(TControlBar *controlBar);
   void            AddControlBar(TControlBar &controlBar);
   void            AddSeparator();
   TControlBarButton *GetClicked() const;
   TControlBarImp *GetControlBarImp() const   { return fControlBarImp; }
   TList          *GetListOfButtons() const   { return fButtons; }
   Int_t           GetNumberOfColumns() const { return fNoroc; }
   Int_t           GetNumberOfRows() const    { return fNoroc; }
   Int_t           GetOrientation() const     { return fOrientation; }
   void            Hide();
   void            SetButtonState(const char *label, Int_t state = 0);
   void            SetFont(const char *fontName);
   void            SetTextColor(const char *colorName);
   void            SetNumberOfColumns(Int_t n) { fNoroc = n; }
   void            SetNumberOfRows(Int_t n) { fNoroc = n; }
   void            SetOrientation(const char *o);
   void            SetOrientation(Int_t o);
   void            SetButtonWidth(UInt_t width);
   void            Show();

   ClassDef(TControlBar,0) //Control bar
};

#endif
