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
/** \class RooStats::HistFactory::Data
 *  \ingroup HistFactory 
*/
 

#include "RooStats/HistFactory/Data.h"


RooStats::HistFactory::Data::Data() : fName("") {
  ;
}

RooStats::HistFactory::Data::Data( std::string HistoName, std::string InputFile, 
				   std::string HistoPath ) :
  fInputFile( InputFile ), fHistoName( HistoName ), fHistoPath( HistoPath ) {;}

TH1* RooStats::HistFactory::Data::GetHisto() {
  return (TH1*) fhData.GetObject();
}

const TH1* RooStats::HistFactory::Data::GetHisto() const {
  return (TH1*) fhData.GetObject();
}


void RooStats::HistFactory::Data::Print( std::ostream& stream ) {


  stream << "\t \t InputFile: " << fInputFile
	 << "\t HistoName: " << fHistoName
	 << "\t HistoPath: " << fHistoPath
	 << "\t HistoAddress: " << GetHisto()
	 << std::endl;

}

void RooStats::HistFactory::Data::writeToFile( std::string OutputFileName, std::string DirName ) {

  TH1* histData = GetHisto();

  if( histData != NULL) {

    histData->Write();
  
    // Set the location of the data
    // in the output measurement
  
    fInputFile = OutputFileName;
    fHistoName = histData->GetName();
    fHistoPath = DirName;

  }

}


void RooStats::HistFactory::Data::PrintXML( std::ostream& xml ) {

  xml << "    <Data HistoName=\"" << GetHistoName() << "\" "
      << "InputFile=\"" << GetInputFile() << "\" "
      << "HistoPath=\"" << GetHistoPath() << "\" "
      << " /> " << std::endl << std::endl;  
  
}
