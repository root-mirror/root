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
/** \class RooStats::HistFactory::Channel
 *  \ingroup HistFactory
  This class encapsulates all information for the statistical interpretation of one experiment.
  It can be combined with other channels (e.g. for the combination of multiple experiments, or
  to constrain nuisance parameters with information obtained in a control region).
  A channel contains one or more samples which describe the contribution from different processes
  to this measurement.
*/



#include "RooStats/HistFactory/Channel.h"
#include "HFMsgService.h"
#include <stdlib.h>

#include "TFile.h"
#include "TKey.h"
#include "TTimeStamp.h"

#include "RooStats/HistFactory/HistFactoryException.h"

using namespace std;

RooStats::HistFactory::Channel::Channel() :
  fName( "" )
{
  // standard constructor
}

RooStats::HistFactory::Channel::Channel(const Channel& other) :
  fName( other.fName ),
  fInputFile( other.fInputFile ),
  fHistoPath( other.fHistoPath ),
  fData( other.fData ),
  fAdditionalData( other.fAdditionalData ),
  fStatErrorConfig( other.fStatErrorConfig ),
  fSamples( other.fSamples )
{ ; }


RooStats::HistFactory::Channel::Channel(std::string ChanName, std::string ChanInputFile) :
  fName( ChanName ), fInputFile( ChanInputFile )
{
  // create channel with given name and input file
}

namespace RooStats{
  namespace HistFactory{
    //BadChannel = Channel();
    Channel BadChannel;
    //    BadChannel.Name = "BadChannel"; // = Channel(); //.Name = "BadChannel";
  }
}


void RooStats::HistFactory::Channel::AddSample( RooStats::HistFactory::Sample sample )
{
  // add fully configured sample to channel
  
  sample.SetChannelName( GetName() );
  fSamples.push_back( sample ); 
}

void RooStats::HistFactory::Channel::Print( std::ostream& stream ) {
  // print information of channel to given stream

  stream << "\t Channel Name: " << fName
	 << "\t InputFile: " << fInputFile
	 << std::endl;

  stream << "\t Data:" << std::endl;
  fData.Print( stream );


  stream << "\t statErrorConfig:" << std::endl;
  fStatErrorConfig.Print( stream );


  if( fSamples.size() != 0 ) {

    stream << "\t Samples: " << std::endl;
    for( unsigned int i = 0; i < fSamples.size(); ++i ) {
      fSamples.at(i).Print( stream );
    }
  }

  
  stream << "\t End of Channel " << fName <<  std::endl;


}  


void RooStats::HistFactory::Channel::PrintXML( std::string directory, std::string prefix ) {

  // Create an XML file for this channel
  cxcoutPHF << "Printing XML Files for channel: " << GetName() << std::endl;
  
  std::string XMLName = prefix + fName + ".xml";
  if( directory != "" ) XMLName = directory + "/" + XMLName;
  
  ofstream xml( XMLName.c_str() );

  // Add the time
  xml << "<!--" << std::endl;
  xml << "This xml file created automatically on: " << std::endl;
  // LM: use TTimeStamp since time_t does not work on Windows
  TTimeStamp t; 
  UInt_t year = 0; 
  UInt_t month = 0; 
  UInt_t day = 0; 
  t.GetDate(true, 0, &year, &month, &day);
  xml << year << '-'
      << month << '-'
      << day
      << std::endl;
  xml << "-->" << std::endl;

  // Add the DOCTYPE
  xml << "<!DOCTYPE Channel  SYSTEM 'HistFactorySchema.dtd'>  " << std::endl << std::endl;

  // Add the Channel
  xml << "  <Channel Name=\"" << fName << "\" InputFile=\"" << fInputFile << "\" >" << std::endl << std::endl;

  fData.PrintXML( xml );
  /*
  xml << "    <Data HistoName=\"" << fData.GetHistoName() << "\" "
      << "InputFile=\"" << fData.GetInputFile() << "\" "
      << "HistoPath=\"" << fData.GetHistoPath() << "\" "
      << " /> " << std::endl << std::endl;  
  */

  fStatErrorConfig.PrintXML( xml );
  /*
  xml << "    <StatErrorConfig RelErrorThreshold=\"" << fStatErrorConfig.GetRelErrorThreshold() << "\" "
      << "ConstraintType=\"" << Constraint::Name( fStatErrorConfig.GetConstraintType() ) << "\" "
      << "/> " << std::endl << std::endl;            
  */

  for( unsigned int i = 0; i < fSamples.size(); ++i ) {
    fSamples.at(i).PrintXML( xml );
    xml << std::endl << std::endl;
  }

  xml << std::endl;
  xml << "  </Channel>  " << std::endl;
  xml.close();

  cxcoutPHF << "Finished printing XML files" << std::endl;

}


void RooStats::HistFactory::Channel::SetData( std::string DataHistoName, std::string DataInputFile, std::string DataHistoPath ) {
  // set data for this channel by specifying the name of the histogram,
  // the external ROOT file and the path to the histogram inside the ROOT file

  fData.SetHistoName( DataHistoName );
  fData.SetInputFile( DataInputFile );
  fData.SetHistoPath( DataHistoPath );

}



void RooStats::HistFactory::Channel::SetData( TH1* hData ) {
  // set data directly to some histogram
  fData.SetHisto( hData ); 
}

void RooStats::HistFactory::Channel::SetData( double val ) {

  // For a NumberCounting measurement only
  // Set the value of data in a particular channel
  // 
  // Internally, this simply creates a 1-bin TH1F for you

  std::string DataHistName = fName + "_data";
  
  // Histogram has 1-bin (hard-coded)
  TH1F* hData = new TH1F( DataHistName.c_str(), DataHistName.c_str(), 1, 0, 1 );
  hData->SetBinContent( 1, val );

  // Set the histogram of the internally held data
  // node of this channel to this newly created histogram
  SetData( hData );

}


void RooStats::HistFactory::Channel::SetStatErrorConfig( double StatRelErrorThreshold, Constraint::Type StatConstraintType ) {

  fStatErrorConfig.SetRelErrorThreshold( StatRelErrorThreshold );
  fStatErrorConfig.SetConstraintType( StatConstraintType );

}

void RooStats::HistFactory::Channel::SetStatErrorConfig( double StatRelErrorThreshold, std::string StatConstraintType ) {

  fStatErrorConfig.SetRelErrorThreshold( StatRelErrorThreshold );
  fStatErrorConfig.SetConstraintType( Constraint::GetType(StatConstraintType) );

}



void RooStats::HistFactory::Channel::CollectHistograms() {

  // Loop through all Samples and Systematics
  // and collect all necessary histograms

  // Handles to open files for collecting histograms
  std::map<std::string,std::unique_ptr<TFile>> fileHandles;

  // Get the Data Histogram:

  if( fData.GetInputFile() != "" ) {
    fData.SetHisto( GetHistogram(fData.GetInputFile(), 
				 fData.GetHistoPath(),
				 fData.GetHistoName(),
				 fileHandles) );
  }

  // Collect any histograms for additional Datasets
  for( unsigned int i=0; i < fAdditionalData.size(); ++i) {
    RooStats::HistFactory::Data& data = fAdditionalData.at(i);
    if( data.GetInputFile() != "" ) {
      data.SetHisto( GetHistogram(data.GetInputFile(), data.GetHistoPath(), data.GetHistoName(), fileHandles) );
    }
  }

  // Get the histograms for the samples:
  for( unsigned int sampItr = 0; sampItr < fSamples.size(); ++sampItr ) {

    RooStats::HistFactory::Sample& sample = fSamples.at( sampItr );


    // Get the nominal histogram:
    cxcoutDHF << "Collecting Nominal Histogram" << std::endl;
    TH1* Nominal =  GetHistogram(sample.GetInputFile(),
				 sample.GetHistoPath(),
				 sample.GetHistoName(),
				 fileHandles);

    sample.SetHisto( Nominal );


    // Get the StatError Histogram (if necessary)
    if( sample.GetStatError().GetUseHisto() ) {
      sample.GetStatError().SetErrorHist( GetHistogram(sample.GetStatError().GetInputFile(),
						       sample.GetStatError().GetHistoPath(),
						       sample.GetStatError().GetHistoName(),
						       fileHandles) );
    }

      
    // Get the HistoSys Variations:
    for( unsigned int histoSysItr = 0; histoSysItr < sample.GetHistoSysList().size(); ++histoSysItr ) {

      RooStats::HistFactory::HistoSys& histoSys = sample.GetHistoSysList().at( histoSysItr );
	
      histoSys.SetHistoLow( GetHistogram(histoSys.GetInputFileLow(), 
					 histoSys.GetHistoPathLow(),
					 histoSys.GetHistoNameLow(),
					 fileHandles) );
	
      histoSys.SetHistoHigh( GetHistogram(histoSys.GetInputFileHigh(),
					  histoSys.GetHistoPathHigh(),
					  histoSys.GetHistoNameHigh(),
					  fileHandles) );
    } // End Loop over HistoSys


      // Get the HistoFactor Variations:
    for( unsigned int histoFactorItr = 0; histoFactorItr < sample.GetHistoFactorList().size(); ++histoFactorItr ) {

      RooStats::HistFactory::HistoFactor& histoFactor = sample.GetHistoFactorList().at( histoFactorItr );

      histoFactor.SetHistoLow( GetHistogram(histoFactor.GetInputFileLow(), 
					    histoFactor.GetHistoPathLow(),
					    histoFactor.GetHistoNameLow(),
					    fileHandles) );
	
      histoFactor.SetHistoHigh( GetHistogram(histoFactor.GetInputFileHigh(),
					     histoFactor.GetHistoPathHigh(),
					     histoFactor.GetHistoNameHigh(),
					     fileHandles) );
    } // End Loop over HistoFactor


      // Get the ShapeSys Variations:
    for( unsigned int shapeSysItr = 0; shapeSysItr < sample.GetShapeSysList().size(); ++shapeSysItr ) {
	
      RooStats::HistFactory::ShapeSys& shapeSys = sample.GetShapeSysList().at( shapeSysItr );

      shapeSys.SetErrorHist( GetHistogram(shapeSys.GetInputFile(), 
					  shapeSys.GetHistoPath(),
					  shapeSys.GetHistoName(),
					  fileHandles) );
    } // End Loop over ShapeSys

    
    // Get any initial shape for a ShapeFactor
    for( unsigned int shapeFactorItr = 0; shapeFactorItr < sample.GetShapeFactorList().size(); ++shapeFactorItr ) {

      RooStats::HistFactory::ShapeFactor& shapeFactor = sample.GetShapeFactorList().at( shapeFactorItr );

      // Check if we need an InitialShape
      if( shapeFactor.HasInitialShape() ) {
	TH1* hist = GetHistogram( shapeFactor.GetInputFile(), shapeFactor.GetHistoPath(), 
				  shapeFactor.GetHistoName(), fileHandles );
	shapeFactor.SetInitialShape( hist );
      }

    } // End Loop over ShapeFactor

  } // End Loop over Samples
}


bool RooStats::HistFactory::Channel::CheckHistograms() const { 

  // Check that all internal histogram pointers
  // are properly configured (ie that they're not NULL)

  try {
  
    if( fData.GetHisto() == NULL && fData.GetInputFile() != "" ) {
      cxcoutEHF << "Error: Data Histogram for channel " << GetName() << " is NULL." << std::endl;
      throw hf_exc();
    }

    // Get the histograms for the samples:
    for( unsigned int sampItr = 0; sampItr < fSamples.size(); ++sampItr ) {

      const RooStats::HistFactory::Sample& sample = fSamples.at( sampItr );

      // Get the nominal histogram:
      if( sample.GetHisto() == NULL ) {
	cxcoutEHF << "Error: Nominal Histogram for sample " << sample.GetName() << " is NULL." << std::endl;
	throw hf_exc();
      } 
      else {

	// Check if any bins are negative
	std::vector<int> NegativeBinNumber;
	std::vector<double> NegativeBinContent;
	const TH1* histNominal = sample.GetHisto();
	for(int ibin=1; ibin<=histNominal->GetNbinsX(); ++ibin) {
	  if(histNominal->GetBinContent(ibin) < 0) {
	    NegativeBinNumber.push_back(ibin);
	    NegativeBinContent.push_back(histNominal->GetBinContent(ibin));
	  }
	}
	if(NegativeBinNumber.size()>0) {
	  cxcoutWHF << "WARNING: Nominal Histogram " << histNominal->GetName() << " for Sample = " << sample.GetName()
		    << " in Channel = " << GetName() << " has negative entries in bin numbers = ";

	  for(unsigned int ibin=0; ibin<NegativeBinNumber.size(); ++ibin) {
	    if(ibin>0) std::cout << " , " ;
	    std::cout << NegativeBinNumber[ibin] << " : " << NegativeBinContent[ibin] ;
	  }
	  std::cout << std::endl;
	}
	
      }

      // Get the StatError Histogram (if necessary)
      if( sample.GetStatError().GetUseHisto() ) {
	if( sample.GetStatError().GetErrorHist() == NULL ) {
	  cxcoutEHF << "Error: Statistical Error Histogram for sample " << sample.GetName() << " is NULL." << std::endl;
	  throw hf_exc();
	}
      }

      
      // Get the HistoSys Variations:
      for( unsigned int histoSysItr = 0; histoSysItr < sample.GetHistoSysList().size(); ++histoSysItr ) {

	const RooStats::HistFactory::HistoSys& histoSys = sample.GetHistoSysList().at( histoSysItr );

	if( histoSys.GetHistoLow() == NULL ) {
	  cxcoutEHF << "Error: HistoSyst Low for Systematic " << histoSys.GetName()
		    << " in sample " << sample.GetName() << " is NULL." << std::endl;
	  throw hf_exc();
	}
	if( histoSys.GetHistoHigh() == NULL ) {
	  cxcoutEHF << "Error: HistoSyst High for Systematic " << histoSys.GetName()
		    << " in sample " << sample.GetName() << " is NULL." << std::endl;
	  throw hf_exc();
	}
	
      } // End Loop over HistoSys


      // Get the HistoFactor Variations:
      for( unsigned int histoFactorItr = 0; histoFactorItr < sample.GetHistoFactorList().size(); ++histoFactorItr ) {

	const RooStats::HistFactory::HistoFactor& histoFactor = sample.GetHistoFactorList().at( histoFactorItr );

	if( histoFactor.GetHistoLow() == NULL ) {
	  cxcoutEHF << "Error: HistoSyst Low for Systematic " << histoFactor.GetName()
		    << " in sample " << sample.GetName() << " is NULL." << std::endl;
	  throw hf_exc();
	}
	if( histoFactor.GetHistoHigh() == NULL ) {
	  cxcoutEHF << "Error: HistoSyst High for Systematic " << histoFactor.GetName()
		    << " in sample " << sample.GetName() << " is NULL." << std::endl;
	  throw hf_exc();
	}

      } // End Loop over HistoFactor


      // Get the ShapeSys Variations:
      for( unsigned int shapeSysItr = 0; shapeSysItr < sample.GetShapeSysList().size(); ++shapeSysItr ) {
	
	const RooStats::HistFactory::ShapeSys& shapeSys = sample.GetShapeSysList().at( shapeSysItr );

	if( shapeSys.GetErrorHist() == NULL ) {
	  cxcoutEHF << "Error: HistoSyst High for Systematic " << shapeSys.GetName()
		    << " in sample " << sample.GetName() << " is NULL." << std::endl;
	  throw hf_exc();
	}

      } // End Loop over ShapeSys

    } // End Loop over Samples

  }
  catch(std::exception& e)
    {
      std::cout << e.what() << std::endl;
      return false;
    }

  return true;




}



/// Open a file and copy a histogram
/// \param InputFile File where the histogram resides.
/// \param HistoPath Path of the histogram in the file.
/// \param HistoName Name of the histogram to retrieve.
/// \param lsof List of open files. Helps to prevent opening and closing a file hundreds of times.
TH1* RooStats::HistFactory::Channel::GetHistogram(std::string InputFile, std::string HistoPath, std::string HistoName, std::map<std::string,std::unique_ptr<TFile>>& lsof) {

  cxcoutPHF << "Getting histogram " << InputFile << ":" << HistoPath << "/" << HistoName << std::endl;

  auto& inFile = lsof[InputFile];
  if (!inFile || !inFile->IsOpen()) {
    inFile.reset( TFile::Open(InputFile.c_str()) );
    if ( !inFile || !inFile->IsOpen() ) {
      cxcoutEHF << "Error: Unable to open input file: " << InputFile << std::endl;
      throw hf_exc();
    }
    cxcoutIHF << "Opened input file: " << InputFile << ": " << std::endl;
  }

  TDirectory* dir = inFile->GetDirectory(HistoPath.c_str());
  if (dir == nullptr) {
    cxcoutEHF << "Histogram path '" << HistoPath
        << "' wasn't found in file '" << InputFile << "'." << std::endl;
    throw hf_exc();
  }

  // Have to read histograms via keys, to ensure that the latest-greatest
  // name cycle is read from file. Otherwise, they might come from memory.
  auto key = dir->GetKey(HistoName.c_str());
  if (key == nullptr) {
    cxcoutEHF << "Histogram '" << HistoName
        << "' wasn't found in file '" << InputFile
        << "' in directory '" << HistoPath << "'." << std::endl;
    throw hf_exc();
  }

  auto hist = dynamic_cast<TH1*>(key->ReadObj());
  if( !hist ) {
    cxcoutEHF << "Histogram '" << HistoName
        << "' wasn't found in file '" << InputFile
        << "' in directory '" << HistoPath << "'." << std::endl;
    throw hf_exc();
  }


  TH1 * ptr = (TH1 *) hist->Clone();

  if(!ptr){
    std::cerr << "Not all necessary info are set to access the input file. Check your config" << std::endl;
    std::cerr << "filename: " << InputFile
	      << "path: " << HistoPath
	      << "obj: " << HistoName << std::endl;
    throw hf_exc();
  }

  ptr->SetDirectory(nullptr);

  
#ifdef DEBUG
  std::cout << "Found Histogram: " << HistoName " at address: " << ptr 
	    << " with integral "   << ptr->Integral() << " and mean " << ptr->GetMean() 
	    << std::endl;
#endif

  // Done
  return ptr;

}
