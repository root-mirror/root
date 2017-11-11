// @(#)root/fitpanel:$Id: LinkDef.h,v 1.0 2003/11/25

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class TFitEditor;
#pragma link C++ class TFitParametersDialog;
#pragma link C++ class TTreeInput;
#pragma link C++ class TAdvancedGraphicsDialog;

#ifdef ROOT7_TFitPanel
#pragma link C++ struct ROOT::Experimental::ComboBoxItem+;
#pragma link C++ class std::vector<ROOT::Experimental::ComboBoxItem>+;
#pragma link C++ struct ROOT::Experimental::TFitPanelModel+;
#pragma link C++ class ROOT::Experimental::TFitPanel+;
#endif

#endif
