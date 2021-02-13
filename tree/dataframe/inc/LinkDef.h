// Author: Danilo Piparo CERN, Enrico Guiraud  4/2018, Vincenzo Eduardo Padulano 06/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                    *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CLING__

#pragma link C++ nestedtypedefs;
#pragma link C++ nestedclasses;

#pragma link C++ class ROOT::RDataFrame-;
#pragma link C++ class ROOT::RDF::RInterface<ROOT::Detail::RDF::RFilterBase>-;
#pragma link C++ class ROOT::RDF::RInterface<ROOT::Detail::RDF::RDefineBase>-;
#pragma link C++ namespace ROOT::Internal::RDF;
#pragma link C++ namespace ROOT::Internal::RDF::GraphDrawing;
#pragma link C++ namespace ROOT::Detail::RDF;
#pragma link C++ namespace ROOT::RDF;
#pragma link C++ class ROOT::Internal::RDF::RActionBase-;
#pragma link C++ class ROOT::Detail::RDF::RFilterBase-;
#pragma link C++ class ROOT::Detail::RDF::RJittedFilter-;
#pragma link C++ class ROOT::Detail::RDF::RDefineBase-;
#pragma link C++ class ROOT::Detail::RDF::RJittedDefine-;
#pragma link C++ class ROOT::Internal::RDF::CountHelper-;
#pragma link C++ class ROOT::Detail::RDF::RRangeBase-;
#pragma link C++ class ROOT::Detail::RDF::RLoopManager-;
#pragma link C++ class ROOT::RDF::TH1DModel-;
#pragma link C++ class ROOT::RDF::TH2DModel-;
#pragma link C++ class ROOT::RDF::TH3DModel-;
#pragma link C++ class ROOT::RDF::TProfile1DModel-;
#pragma link C++ class ROOT::RDF::TProfile2DModel-;
#pragma link C++ class ROOT::Internal::RDF::RIgnoreErrorLevelRAII-;
#pragma link C++ class ROOT::Internal::RDF::FillHelper-;
#pragma link C++ class ROOT::RDF::RTrivialDS-;
#pragma link C++ class ROOT::Internal::RDF::RRootDS-;
#pragma link C++ class ROOT::RDF::RCsvDS-;
#pragma link C++ class ROOT::Internal::RDF::MeanHelper-;
#pragma link C++ class ROOT::Internal::RDF::RBookedDefines-;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValueBase+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<int>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<unsigned int>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<float>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<double>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<Long64_t>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<ULong64_t>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<Double_t>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<TH1D>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<TH2D>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<TH3D>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<TGraph>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<TStatistic>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<TProfile>+;
#pragma link C++ class ROOT::Detail::RDF::RMergeableValue<TProfile2D>+;

#endif


