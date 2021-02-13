// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, George Lewis
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////



#include <string>
#include <exception>
#include <vector>

#include <TROOT.h>

//void topDriver(string input); // in MakeModelAndMeasurements
//void fastDriver(string input); // in MakeModelAndMeasurementsFast

//#include "RooStats/HistFactory/MakeModelAndMeasurements.h"
#include "RooStats/HistFactory/ConfigParser.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "HFMsgService.h"
#include "hist2workspaceCommandLineOptionsHelp.h"

//_____________________________batch only_____________________
#ifndef __CINT__

namespace RooStats {
  namespace HistFactory {
    void fastDriver(std::string input) {

      // Create the initial list of measurements and channels
      std::vector< HistFactory::Measurement > measurement_list;
      // std::vector< HistFactory::Channel >     channel_list;

      // Fill them using the XML parser
      HistFactory::ConfigParser xmlParser;
      measurement_list = xmlParser.GetMeasurementsFromXML( input );

      // At this point, we have all the information we need
      // from the xml files.

      // We will make the measurements 1-by-1
      for(unsigned int i = 0; i < measurement_list.size(); ++i) {
	HistFactory::Measurement measurement = measurement_list.at(i);
	measurement.CollectHistograms();
	MakeModelAndMeasurementFast( measurement );
      }

      return;

    }
  } // namespace RooStats
} // namespace HistFactory

/**
 * \ingroup HistFactory
 * main function of the hist2workspace executable.
 * It creates RooFit models from an xml config and files with histograms.
 * See MakeModelAndMeasurementFast(), for further instructions.
 * \param[in] -h Help
 * \param[in] -standard_form Standard xml model definitions. See MakeModelAndMeasurementFast()
 * \param[in] -number_counting_form Deprecated
 * \param[in] -v Switch HistFactory message stream to INFO level.
 * \param[in] -vv Switch HistFactory message stream to DEBUG level.
 */
int main(int argc, char** argv) {

  if( !(argc > 1) ) {
    std::cerr << "need input file" << std::endl;
    exit(1);
  }

  //Switch off ROOT histogram memory management
  gROOT->SetMustClean(false);
  TDirectory::AddDirectory(false);
  cxcoutIHF << "hist2workspace is less verbose now. Use -v and -vv for more details." << std::endl;
  RooMsgService::instance().getStream(1).minLevel = RooFit::PROGRESS;
  RooMsgService::instance().getStream(2).minLevel = RooFit::PROGRESS;
  RooMsgService::instance().getStream(2).addTopic(RooFit::HistFactory);

  std::string driverArg;

  for (int i=1; i < argc; ++i) {
    std::string input = argv[i];

    if (input == "-h" || input == "--help"){
      fprintf(stderr, kCommandLineOptionsHelp);
      return 0;
    }

    if (input == "-v") {
      RooMsgService::instance().getStream(1).minLevel = RooFit::INFO;
      RooMsgService::instance().getStream(2).minLevel = RooFit::INFO;
      continue;
    }

    if (input == "-vv") {
      RooMsgService::instance().getStream(1).minLevel = RooFit::INFO;
      RooMsgService::instance().getStream(2).minLevel = RooFit::DEBUG;
      continue;
    }

    if (input == "-number_counting_form") {
      std::cout << "ERROR: 'number_counting_form' is now deprecated." << std::endl;
      return 255;
    }

    if(input == "-standard_form") {
      driverArg = argv[++i];
      continue;
    }

    driverArg = argv[i];
  }

  try {
    RooStats::HistFactory::fastDriver(driverArg);
  }
  catch(const std::string &str) {
    std::cerr << "hist2workspace - Caught exception: " << str << std::endl ;
    return 1;
  }
  catch( const std::exception& e ) {
    std::cerr << "hist2workspace - Caught Exception: " << e.what() << std::endl;
    return 1;
  }
  catch(...) {
    std::cerr << "hist2workspace - Caught Exception" << std::endl;
    return 1;
  }


  return 0;
}

#endif
