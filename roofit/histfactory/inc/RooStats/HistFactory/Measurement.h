// @(#)root/roostats:$Id$
// Author: George Lewis, Kyle Cranmer
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef HISTFACTORY_MEASUREMENT_H
#define HISTFACTORY_MEASUREMENT_H

#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <vector>

#include "TNamed.h"
#include "TFile.h"

#include "PreprocessFunction.h"
#include "RooStats/HistFactory/Channel.h"
#include "RooStats/HistFactory/Asimov.h"

namespace RooStats{
namespace HistFactory {

class Measurement : public TNamed {

public:

  Measurement();
  ///  Measurement( const Measurement& other ); // Copy
  Measurement(const char* Name, const char* Title="");

  ///  set output prefix
  void SetOutputFilePrefix( const std::string& prefix ) { fOutputFilePrefix = prefix; }
  /// retrieve prefix for output files
  std::string GetOutputFilePrefix() { return fOutputFilePrefix; }

  /// insert PoI at beginning of vector of PoIs
  void SetPOI( const std::string& POI ) { fPOI.insert( fPOI.begin(), POI ); }
  /// append parameter to vector of PoIs
  void AddPOI( const std::string& POI ) { fPOI.push_back(POI); }
  /// get name of PoI at given index
  std::string GetPOI(unsigned int i=0) { return fPOI.at(i); }
  /// get vector of PoI names
  std::vector<std::string>& GetPOIList() { return fPOI; }


  /// Add a parameter to be set as constant
  /// (Similar to ParamSetting method below)
  void AddConstantParam( const std::string& param );
  /// empty vector of constant parameters
  void ClearConstantParams() { fConstantParams.clear(); }
  /// get vector of all constant parameters
  std::vector< std::string >& GetConstantParams() { return fConstantParams; }

  /// Set a parameter to a specific value
  /// (And optionally fix it)
  void SetParamValue( const std::string& param, double value);
  /// get map: parameter name <--> parameter value
  std::map<std::string, double>& GetParamValues() { return fParamValues; }
  /// clear map of parameter values
  void ClearParamValues() { fParamValues.clear(); }

  void AddPreprocessFunction( std::string name, std::string expression, std::string dependencies );
  /// add a preprocess function object
  void AddFunctionObject( const RooStats::HistFactory::PreprocessFunction function) { fFunctionObjects.push_back( function ); }
  void SetFunctionObjects( std::vector< RooStats::HistFactory::PreprocessFunction > objects ) { fFunctionObjects = objects; }
  /// get vector of defined function objects
  std::vector< RooStats::HistFactory::PreprocessFunction >& GetFunctionObjects() { return fFunctionObjects; }
  const std::vector< RooStats::HistFactory::PreprocessFunction >& GetFunctionObjects() const { return fFunctionObjects; }  
  std::vector< std::string > GetPreprocessFunctions() const;

  /// get vector of defined Asimov Datasets
  std::vector< RooStats::HistFactory::Asimov >& GetAsimovDatasets() { return fAsimovDatasets; }
  /// add an Asimov Dataset
  void AddAsimovDataset( RooStats::HistFactory::Asimov dataset ) { fAsimovDatasets.push_back(dataset); }

  /// set integrated luminosity used to normalise histograms (if NormalizeByTheory is true for this sample)
  void SetLumi(double Lumi ) { fLumi = Lumi; }
  /// set relative uncertainty on luminosity
  void SetLumiRelErr( double RelErr ) { fLumiRelErr = RelErr; }
  /// retrieve integrated luminosity
  double GetLumi() { return fLumi; }
  /// retrieve relative uncertainty on luminosity
  double GetLumiRelErr() { return fLumiRelErr; }
  
  void SetBinLow( int BinLow ) { fBinLow = BinLow; }
  void SetBinHigh ( int BinHigh ) { fBinHigh = BinHigh; }
  int GetBinLow() { return fBinLow; }
  int GetBinHigh() { return fBinHigh; } 

  /// do not produce any plots or tables, just save the model
  void SetExportOnly( bool ExportOnly ) { fExportOnly = ExportOnly; }
  bool GetExportOnly() { return fExportOnly; }

  void PrintTree( std::ostream& = std::cout ); /// Print to a stream
  void PrintXML( std::string Directory="", std::string NewOutputPrefix="" );

  std::vector< RooStats::HistFactory::Channel >& GetChannels() { return fChannels; }
  const std::vector< RooStats::HistFactory::Channel >& GetChannels() const { return fChannels; } 
  RooStats::HistFactory::Channel& GetChannel( std::string );
  /// add a completely configured channel
  void AddChannel( RooStats::HistFactory::Channel chan ) { fChannels.push_back( chan ); }

  bool HasChannel( std::string );
  void writeToFile( TFile* file );

  void CollectHistograms();


  void AddGammaSyst(std::string syst, double uncert);
  void AddLogNormSyst(std::string syst, double uncert);
  void AddUniformSyst(std::string syst);
  void AddNoSyst(std::string syst);

  std::map< std::string, double >& GetGammaSyst() { return fGammaSyst; }
  std::map< std::string, double >& GetUniformSyst() { return fUniformSyst; }
  std::map< std::string, double >& GetLogNormSyst() { return fLogNormSyst; }
  std::map< std::string, double >& GetNoSyst() { return fNoSyst; }

  std::string GetInterpolationScheme() { return fInterpolationScheme; }
  
private:

  /// Configurables of this measurement
  std::string fOutputFilePrefix;
  std::vector<std::string> fPOI;
  double fLumi;
  double fLumiRelErr;
  int fBinLow;
  int fBinHigh;
  bool fExportOnly;
  std::string fInterpolationScheme;

  /// Channels that make up this measurement
  std::vector< RooStats::HistFactory::Channel > fChannels;

  /// List of Parameters to be set constant
  std::vector< std::string > fConstantParams;

  /// Map of parameter names to inital values to be set
  std::map< std::string, double > fParamValues;

  /// List of Preprocess Function objects
  std::vector< RooStats::HistFactory::PreprocessFunction > fFunctionObjects;

  /// List of Asimov datasets to generate
  std::vector< RooStats::HistFactory::Asimov > fAsimovDatasets;

  /// List of Alternate constraint terms
  std::map< std::string, double > fGammaSyst;
  std::map< std::string, double > fUniformSyst;
  std::map< std::string, double > fLogNormSyst;
  std::map< std::string, double > fNoSyst;

  std::string GetDirPath( TDirectory* dir );

  ClassDef(RooStats::HistFactory::Measurement, 3);

};

} // namespace HistFactory
} // namespace RooStats

#endif
