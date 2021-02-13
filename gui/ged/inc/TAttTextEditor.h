// @(#)root/ged:$Id$
// Author: Ilka  Antcheva 11/05/04

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TAttTextEditor
#define ROOT_TAttTextEditor

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TAttTextEditor                                                      //
//                                                                      //
//  Implements GUI for editing text attributes.                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGedFrame.h"

#include "TGSlider.h"

class TGComboBox;
class TGFontTypeComboBox;
class TGColorSelect;
class TAttText;
class TGNumberEntryField;

class TAttTextEditor : public TGedFrame {

protected:
   TAttText            *fAttText;         // text attribute object
   TGFontTypeComboBox  *fTypeCombo;       // font style combo box
   TGComboBox          *fSizeCombo;       // font size combo box
   TGComboBox          *fAlignCombo;      // font aligh combo box
   TGColorSelect       *fColorSelect;     // color selection widget
   TGHSlider           *fAlpha;           // fill opacity
   TGNumberEntryField  *fAlphaField;

   void             ConnectSignals2Slots();

   static  TGComboBox *BuildFontSizeComboBox(TGFrame *parent, Int_t id);
   static  TGComboBox *BuildTextAlignComboBox(TGFrame *parent, Int_t id);

public:
   TAttTextEditor(const TGWindow *p = 0,
                  Int_t width = 140, Int_t height = 30,
                  UInt_t options = kChildFrame,
                  Pixel_t back = GetDefaultFrameBackground());
   virtual ~TAttTextEditor();

   virtual void     SetModel(TObject* obj);
   virtual Bool_t   ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
   virtual void     DoTextAlphaColor(ULong_t p);
   virtual void     DoAlpha();
   virtual void     DoAlphaField();
   virtual void     DoLiveAlpha(Int_t a);
   virtual void     GetCurAlpha();
   virtual void     DoTextColor(Pixel_t color);

   ClassDef(TAttTextEditor,0)  //GUI for editing text attributes
};

#endif
