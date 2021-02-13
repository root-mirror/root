// @(#)root/roostats:$Id$
// Author: George Lewis, Kyle Cranmer
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef HISTFACTORY_CHANNEL_H
#define HISTFACTORY_CHANNEL_H

#include "RooStats/HistFactory/Data.h"
#include "RooStats/HistFactory/Sample.h"
#include "RooStats/HistFactory/Systematics.h"

#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

class TFile;

namespace RooStats{
namespace HistFactory {

class Channel  {


public:
  friend class Measurement;

  Channel();
  Channel(const Channel& other);
  Channel(std::string Name, std::string InputFile="");

  /// set name of channel
  void SetName( const std::string& Name ) { fName = Name; }
  /// get name of channel
  std::string GetName() const { return fName; }
  /// set name of input file containing histograms
  void SetInputFile( const std::string& file ) { fInputFile = file; }
  /// get name of input file
  std::string GetInputFile() const { return fInputFile; }
  /// set path for histograms in input file
  void SetHistoPath( const std::string& file ) { fHistoPath = file; }
  /// get path to histograms in input file
  std::string GetHistoPath() const { return fHistoPath; }

  /// set data object
  void SetData( const RooStats::HistFactory::Data& data ) { fData = data; }
  void SetData( std::string HistoName, std::string InputFile, std::string HistoPath="" );
  void SetData( double Val );
  void SetData( TH1* hData );
  /// get data object
  RooStats::HistFactory::Data& GetData() { return fData; }

  /// add additional data object
  void AddAdditionalData( const RooStats::HistFactory::Data& data ) { fAdditionalData.push_back(data); }
  /// retrieve vector of additional data objects
  std::vector<RooStats::HistFactory::Data>& GetAdditionalData() { return fAdditionalData; }

  void SetStatErrorConfig( double RelErrorThreshold, Constraint::Type ConstraintType );
  void SetStatErrorConfig( double RelErrorThreshold, std::string ConstraintType );
  /// define treatment of statistical uncertainties
  void SetStatErrorConfig( RooStats::HistFactory::StatErrorConfig Config ) { fStatErrorConfig = Config; }
  /// get information about threshold for statistical uncertainties and constraint term
  HistFactory::StatErrorConfig& GetStatErrorConfig() { return fStatErrorConfig; }
  const HistFactory::StatErrorConfig& GetStatErrorConfig() const { return fStatErrorConfig; }  

  void AddSample( RooStats::HistFactory::Sample sample );
  /// get vector of samples for this channel
  std::vector< RooStats::HistFactory::Sample >& GetSamples() { return fSamples; }
  const std::vector< RooStats::HistFactory::Sample >& GetSamples() const { return fSamples; }

  void Print(std::ostream& = std::cout);  
  void PrintXML( std::string Directory, std::string Prefix="" );
  
  void CollectHistograms();
  bool CheckHistograms() const;

protected:
  
  std::string fName;
  std::string fInputFile;
  std::string fHistoPath;

  HistFactory::Data fData;

  /// One can add additional datasets
  /// These are simply added to the xml under a different name
  std::vector<RooStats::HistFactory::Data> fAdditionalData;

  HistFactory::StatErrorConfig fStatErrorConfig;

  std::vector< RooStats::HistFactory::Sample > fSamples;

  TH1* GetHistogram( std::string InputFile, std::string HistoPath, std::string HistoName, std::map<std::string, std::unique_ptr<TFile>>& lsof);
};

  extern Channel BadChannel;
 
} // namespace HistFactory
} // namespace RooStats

#endif
