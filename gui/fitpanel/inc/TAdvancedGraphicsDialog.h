// @(#)root/fitpanel:$Id$
// Author: David Gonzalez Maline 11/12/2008

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT__TAdvancedGraphicsDialog__
#define ROOT__TAdvancedGraphicsDialog__

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAdvancedGraphicsDialog                                              //
//                                                                      //
// Allows to create advanced graphics from the last fit made in the     //
// fitpanel. This includes the scan graphics, the contour and the       //
// confidence levels.                                                   //
//////////////////////////////////////////////////////////////////////////

#include "TGFrame.h"
#include "TGButton.h"
#include "TGComboBox.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGTab.h"
#include "TGColorSelect.h"

#include "TBackCompFitter.h"
#include "TF1.h"

enum EAdvanceGraphicsDialog {
   kAGD_TMETHOD,  kAGD_CONTOURMETHOD, kAGD_SCANMETHOD,
   kAGD_CONTPAR1, kAGD_CONTPAR2,      kAGD_CONTERR,
   kAGD_CONTOVER, kAGD_CONTCOLOR,
   kAGD_BDRAW, kAGD_BCLOSE,
   kAGD_SCANPAR, kAGD_SCANMIN, kAGD_SCANMAX,

   kAGD_PARCOUNTER = 1000
};

class TAdvancedGraphicsDialog : public TGTransientFrame {

private:
   TGVerticalFrame  *fMainFrame;     // Main Vertical Frame
   TGTab            *fTab;           // Tab containing the available methods

   TGVerticalFrame  *fContourFrame;  // Contour Frame
   TGNumberEntry    *fContourPoints; // Number of points for the graph
   TGComboBox       *fContourPar1;   // Parameter 1 for Contour
   TGComboBox       *fContourPar2;   // Parameter 2 for Contour
   TGNumberEntry    *fContourError;  // Error Level for Contour
   TGCheckButton    *fContourOver;   // Superimpose the graphics
   TGColorSelect    *fContourColor;  // Color for the graph

   TGVerticalFrame  *fScanFrame;     // Scan Frame
   TGNumberEntry    *fScanPoints;    // Number of points for the graph
   TGComboBox       *fScanPar;       // Parameter for Scan
   TGNumberEntry    *fScanMin;       // Min Value for Contour
   TGNumberEntry    *fScanMax;       // Max Value for Contour

   TGVerticalFrame  *fConfFrame;     // Confidence Intervals Frame
   TGNumberEntry    *fConfLevel;     // Confidence Level
   TGColorSelect    *fConfColor;     // Color for the graph

   TGTextButton     *fDraw;          // ok button
   TGTextButton     *fClose;         // cancel button

   TBackCompFitter  *fFitter;        // Fitter.

   void CreateContourFrame();
   void CreateScanFrame();
   void CreateConfFrame();
   void AddParameters(TGComboBox*);

   void DrawContour();
   void DrawScan();
   void DrawConfidenceLevels();

   void ConnectSlots();

   TAdvancedGraphicsDialog(const TAdvancedGraphicsDialog&);  // Not implemented
   TAdvancedGraphicsDialog &operator= (const TAdvancedGraphicsDialog&); // Not implemented

public:
   TAdvancedGraphicsDialog(const TGWindow *p, const TGWindow *main);
   ~TAdvancedGraphicsDialog();

   void DoDraw();
   void DoChangedScanPar(Int_t selected);

   ClassDef(TAdvancedGraphicsDialog, 0)  // Simple input dialog
};

#endif

