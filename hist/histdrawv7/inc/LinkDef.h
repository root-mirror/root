/* @(#)root/hist:$Id$ */

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// need to correctly generate dictionary for display items, used in v7 histpainter,
// but currently histpainter does not creates dictionary at all
// #pragma extra_include "ROOT/RDisplayItem.hxx";

#pragma link C++ class ROOT::Experimental::Internal::RIOShared<ROOT::Experimental::Detail::RHistImplPrecisionAgnosticBase<1>>+;
#pragma link C++ class ROOT::Experimental::Internal::RIOShared<ROOT::Experimental::Detail::RHistImplPrecisionAgnosticBase<2>>+;
#pragma link C++ class ROOT::Experimental::Internal::RIOShared<ROOT::Experimental::Detail::RHistImplPrecisionAgnosticBase<3>>+;

#pragma link C++ class ROOT::Experimental::RHistDrawableBase+;
#pragma link C++ class ROOT::Experimental::RHistDrawableBase::RRequest+;
#pragma link C++ class ROOT::Experimental::RHistDrawableBase::RReply+;
#pragma link C++ class ROOT::Experimental::RHistDrawable<1>+;
#pragma link C++ class ROOT::Experimental::RHistDrawable<2>+;
#pragma link C++ class ROOT::Experimental::RHistDrawable<3>+;
#pragma link C++ class ROOT::Experimental::RHist1Drawable+;
#pragma link C++ class ROOT::Experimental::RHist2Drawable+;
#pragma link C++ class ROOT::Experimental::RHist3Drawable+;


#pragma link C++ class ROOT::Experimental::RHistDisplayItem+;

#pragma link C++ class ROOT::Experimental::RDisplayHistStat+;

#pragma link C++ class ROOT::Experimental::RHistStatBoxBase+;
#pragma link C++ class ROOT::Experimental::RHistStatBoxBase::RRequest+;
#pragma link C++ class ROOT::Experimental::RHistStatBoxBase::RReply+;

#pragma link C++ class ROOT::Experimental::RHistStatBox<1>+;
#pragma link C++ class ROOT::Experimental::RHistStatBox<2>+;
#pragma link C++ class ROOT::Experimental::RHistStatBox<3>+;

#pragma link C++ class ROOT::Experimental::RHist1StatBox+;
#pragma link C++ class ROOT::Experimental::RHist2StatBox+;
#pragma link C++ class ROOT::Experimental::RHist3StatBox+;

#endif
