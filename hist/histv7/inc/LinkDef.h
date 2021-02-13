/* @(#)root/hist:$Id$ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class ROOT::Experimental::RH1F+;
#pragma link C++ class ROOT::Experimental::RH1D+;
#pragma link C++ class ROOT::Experimental::RH2F+;
#pragma link C++ class ROOT::Experimental::RH2D+;
#pragma link C++ class ROOT::Experimental::RH3F+;
#pragma link C++ class ROOT::Experimental::RH3D+;
// #pragma link C++ class ROOT::Experimental::Detail::RHistImpl<ROOT::Experimental::Detail::RHistData<1,double,vector<double>,ROOT::Experimental::RHistStatContent,ROOT::Experimental::RHistStatUncertainty>,ROOT::Experimental::RAxisEquidistant>+;
// #pragma link C++ class ROOT::Experimental::Detail::RHistImpl<ROOT::Experimental::Detail::RHistData<2,double,vector<double>, ROOT::Experimental::RHistStatContent, ROOT::Experimental::RHistStatUncertainty>, ROOT::Experimental::RAxisEquidistant, ROOT::Experimental::RAxisIrregular>+;
#pragma link C++ class ROOT::Experimental::Detail::RHistImplBase<ROOT::Experimental::Detail::RHistData<1,double,vector<double>,ROOT::Experimental::RHistStatContent,ROOT::Experimental::RHistStatUncertainty>>+;
#pragma link C++ class ROOT::Experimental::Detail::RHistImplBase<ROOT::Experimental::Detail::RHistData<2,double,vector<double>,ROOT::Experimental::RHistStatContent,ROOT::Experimental::RHistStatUncertainty>>+;
#pragma link C++ class ROOT::Experimental::Detail::RHistImplBase<ROOT::Experimental::Detail::RHistData<3,double,vector<double>,ROOT::Experimental::RHistStatContent,ROOT::Experimental::RHistStatUncertainty>>+;
#pragma link C++ class ROOT::Experimental::Detail::RHistImplPrecisionAgnosticBase<1>+;
#pragma link C++ class ROOT::Experimental::Detail::RHistImplPrecisionAgnosticBase<2>+;
#pragma link C++ class ROOT::Experimental::Detail::RHistImplPrecisionAgnosticBase<3>+;
#pragma link C++ class ROOT::Experimental::RHistStatContent<1,double>+;
#pragma link C++ class ROOT::Experimental::RHistStatContent<2,double>+;
#pragma link C++ class ROOT::Experimental::RHistStatContent<3,double>+;
#pragma link C++ class ROOT::Experimental::RHistStatUncertainty<1,double>+;
#pragma link C++ class ROOT::Experimental::RHistStatUncertainty<2,double>+;
#pragma link C++ class ROOT::Experimental::RHistStatUncertainty<3,double>+;
#pragma link C++ class ROOT::Experimental::Detail::RHistData<1,double,vector<double>,ROOT::Experimental::RHistStatContent,ROOT::Experimental::RHistStatUncertainty>+;
#pragma link C++ class ROOT::Experimental::Detail::RHistData<2,double,vector<double>,ROOT::Experimental::RHistStatContent,ROOT::Experimental::RHistStatUncertainty>+;
#pragma link C++ class ROOT::Experimental::Detail::RHistData<3,double,vector<double>,ROOT::Experimental::RHistStatContent,ROOT::Experimental::RHistStatUncertainty>+;
#pragma link C++ class ROOT::Experimental::RAxisBase+;
#pragma link C++ class ROOT::Experimental::RAxisEquidistant+;
#pragma link C++ class ROOT::Experimental::RAxisIrregular+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisEquidistant>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisIrregular>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisEquidistant,ROOT::Experimental::RAxisIrregular>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisIrregular,ROOT::Experimental::RAxisEquidistant>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisEquidistant,ROOT::Experimental::RAxisEquidistant>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisIrregular,ROOT::Experimental::RAxisIrregular>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisEquidistant,ROOT::Experimental::RAxisEquidistant,ROOT::Experimental::RAxisEquidistant>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisEquidistant,ROOT::Experimental::RAxisEquidistant,ROOT::Experimental::RAxisIrregular>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisEquidistant,ROOT::Experimental::RAxisIrregular,ROOT::Experimental::RAxisEquidistant>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisEquidistant,ROOT::Experimental::RAxisIrregular,ROOT::Experimental::RAxisIrregular>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisIrregular,ROOT::Experimental::RAxisEquidistant,ROOT::Experimental::RAxisEquidistant>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisIrregular,ROOT::Experimental::RAxisEquidistant,ROOT::Experimental::RAxisIrregular>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisIrregular,ROOT::Experimental::RAxisIrregular,ROOT::Experimental::RAxisEquidistant>+;
#pragma link C++ class tuple<ROOT::Experimental::RAxisIrregular,ROOT::Experimental::RAxisIrregular,ROOT::Experimental::RAxisIrregular>+;
#endif
