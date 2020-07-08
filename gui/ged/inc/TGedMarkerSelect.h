// @(#)root/ged:$Id$
// Author: Marek Biskup, Ilka Antcheva   24/07/03

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGedMarkerSelect
#define ROOT_TGedMarkerSelect

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGedMarkerPopup and TGedMarkerSelect.                                //
//                                                                      //
// The TGedMarkerPopup is a popup containing all diferent styles of     //
// markers.                                                             //
//                                                                      //
// The TGedMarkerSelect widget is a button with marker drawn inside     //
// and a little down arrow. When clicked the TGMarkerPopup.             //
//                                                                      //
// Selecting a marker in this widget will generate the event:           //
// kC_MARKERSEL, kMAR_SELCHANGED, widget id, pixel.                     //
//                                                                      //
// and the signal:                                                      //
// MarkerSelected(Style_t marker)                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGedPatternSelect.h"

class TGToolTip;
class TGPicture;


class TGedMarkerPopup : public TGedPopup {

protected:
   Style_t  fCurrentStyle;     //currently selected style

public:
   TGedMarkerPopup(const TGWindow *p, const TGWindow *m, Style_t markerStyle);
   ~TGedMarkerPopup() override;

   Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) override;

   ClassDef(TGedMarkerPopup,0)  //marker select popup
};


class TGedMarkerSelect : public TGedSelect {

protected:
   Style_t          fMarkerStyle;   //marker style
   const TGPicture *fPicture;       //image used for popup window

   void     DoRedraw() override;

public:
   TGedMarkerSelect(const TGWindow *p, Style_t markerStyle, Int_t id);
   ~TGedMarkerSelect() override { if(fPicture) gClient->FreePicture(fPicture);}

   Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) override;
   Style_t        GetMarkerStyle() const { return fMarkerStyle; }
   void           SetMarkerStyle(Style_t pattern);
   virtual void   MarkerSelected(Style_t marker = 0) { Emit("MarkerSelected(Style_t)", marker ? marker : GetMarkerStyle()); }  // *SIGNAL*
   void   SavePrimitive(std::ostream &out, Option_t * = "") override;
   TGDimension GetDefaultSize() const override { return TGDimension(38, 21); }

   ClassDef(TGedMarkerSelect,0)  // Marker selection button
};

#endif
