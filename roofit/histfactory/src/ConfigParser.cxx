// @(#)root/roostats:$Id:  cranmer $
// Author: Kyle Cranmer, Akira Shibata
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/** \class RooStats::HistFactory::ConfigParser
 *  \ingroup HistFactory 
 *  TODO Add documentation.
*/

#include "TDOMParser.h"

#include "RooStats/HistFactory/ConfigParser.h"
#include "RooStats/HistFactory/HistFactoryException.h"
#include "RooStats/HistFactory/Measurement.h"

#include "Helper.h"
#include "HFMsgService.h"

using namespace RooStats;
using namespace HistFactory;

using namespace std;

std::vector< RooStats::HistFactory::Measurement > ConfigParser::GetMeasurementsFromXML( string input ) {
 
  // Open an input "Driver" XML file (input),
  // Parse that file and its channel files
  // and return a vector filled with 
  // the listed measurements


  // Create the list of measurements
  // (This is what will be returned)
  std::vector< HistFactory::Measurement > measurement_list;

  try {

    // Open the Driver XML File
    TDOMParser xmlparser;
    Int_t parseError = xmlparser.ParseFile( input.c_str() );
    if( parseError ) { 
      std::cerr << "Loading of xml document \"" << input
		<< "\" failed" << std::endl;
      throw hf_exc();
    } 


    // Read the Driver XML File
    cxcoutIHF << "reading input : " << input << endl;
    TXMLDocument* xmldoc = xmlparser.GetXMLDocument();
    TXMLNode* rootNode = xmldoc->GetRootNode();


    // Check that it is the proper DOCTYPE
    if( rootNode->GetNodeName() != TString( "Combination" ) ){
      cxcoutEHF << "Error: Driver DOCTYPE not equal to 'Combination'" << std::endl;
      throw hf_exc();
    }

    // Loop over the Combination's attributes
    std::string OutputFilePrefix;

    TListIter attribIt = rootNode->GetAttributes();
    TXMLAttr* curAttr = 0;
    while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

      // Get the Name, Val of this node
      TString attrName    = curAttr->GetName();
      std::string attrVal = curAttr->GetValue();

      if( attrName == TString( "" ) ) {
	cxcoutEHF << " Error: Attribute for 'Combination' with no name found" << std::endl;
	throw hf_exc();
      }

      else if( attrName == TString( "OutputFilePrefix" ) ) {
	OutputFilePrefix = string(curAttr->GetValue());
	cxcoutIHF << "output file prefix is : " << OutputFilePrefix << endl;
      }

      /*
	else if( attrName == TString( "InputFile" ) ) {
        channel.InputFile = attrVal ;
	}
      */

      else {
	cxcoutEHF << " Error: Unknown attribute for 'Combination' encountered: "
		  << attrName << std::endl;
	throw hf_exc();
      }

      // node = node->GetNextNode();

    }

    TXMLNode* node = NULL;

    // Get the list of channel XML files to combine
    // Do this first so we can quickly exit
    // if no channels are found
    std::vector< std::string > xml_channel_files;
    node = rootNode->GetChildren();
    while( node != 0 ) {
      if( node->GetNodeName() == TString( "Input" ) ) {
	if( node->GetText() == NULL ) {
	  cxcoutEHF << "Error: node: " << node->GetName()
		    << " has no text." << std::endl;
	  throw hf_exc();
	}
	xml_channel_files.push_back(node->GetText());
      }
      node = node->GetNextNode();
    }

    // If no channel xml files are found, exit
    if(xml_channel_files.empty()){
      cerr << "no input channels found" << endl;
      throw hf_exc();
      //return measurement_list;
    }
    else {
      std::ostringstream msg;
      msg << "Found Channels: ";
      for( unsigned int i=0; i < xml_channel_files.size(); ++i )   msg << " " << xml_channel_files.at(i);
      cxcoutIHF << msg.str();
    }

    // Get the list of functions
    // These apply to all measurements, so we
    // first create the list of preprocess functions 
    // (before we create the list of measurements)
    // and then we add them to all measurements

    // For now, we create this list twice
    // simply for compatability
    // std::vector< std::string > preprocessFunctions;
    std::vector< RooStats::HistFactory::PreprocessFunction > functionObjects;

    node = rootNode->GetChildren();
    while( node != 0 ) {
      if( node->GetNodeName() == TString( "Function" ) ) {
      
	// For now, add both the objects itself and
	// it's command string (for easy compatability)
	RooStats::HistFactory::PreprocessFunction Func = ParseFunctionConfig( node );
	// preprocessFunctions.push_back( Func.GetCommand() ); 
	functionObjects.push_back( Func );
      }
      node = node->GetNextNode();
    }

    std::cout << std::endl;


    // Fill the list of measurements
    node = rootNode->GetChildren();
    while( node != 0 ) {

      if( node->GetNodeName() == TString( "" ) ) {
	cxcoutEHF << "Error: Node found in Measurement Driver XML with no name" << std::endl;
	throw hf_exc();
      }

      else if( node->GetNodeName() == TString( "Measurement" ) ) {
	HistFactory::Measurement measurement = CreateMeasurementFromDriverNode( node );
	// Set the prefix (obtained above)
	measurement.SetOutputFilePrefix( OutputFilePrefix );
	measurement_list.push_back( measurement );
      }

      else if( node->GetNodeName() == TString( "Function" ) ) {
	// Already processed these (directly above)
	;
      }

      else if( node->GetNodeName() == TString( "Input" ) ) {
	// Already processed these (directly above)
	;
      }

      else if( IsAcceptableNode( node ) ) { ; }
    
      else {
	cxcoutEHF << "Error: Unknown node found in Measurement Driver XML: "
		  << node->GetNodeName() << std::endl;
	throw hf_exc();
      }

      node = node->GetNextNode();
    }

    cxcoutIHF << "Done Processing Measurements" << std::endl;

    if( measurement_list.size() == 0 ) {
      cxcoutEHF << "Error: No Measurements found in XML Driver File" << std::endl;
      throw hf_exc();
    }
    else {
      std::ostringstream msg;
      msg << "Found Measurements: ";
      for( unsigned int i=0; i < measurement_list.size(); ++i )   msg << " " << measurement_list.at(i).GetName();
      msg << std::endl;
      cxcoutIHF << msg.str();
    }

    // Add the preprocessed functions to each measurement
    // for( unsigned int i = 0; i < measurement_list.size(); ++i) {
    //   measurement_list.at(i).SetPreprocessFunctions( preprocessFunctions );
    // }
    // Add the preprocessed functions to each measurement
    for( unsigned int i = 0; i < measurement_list.size(); ++i) {
      measurement_list.at(i).SetFunctionObjects( functionObjects );
    }

    // Create an instance of the class
    // that collects histograms
    //HistCollector collector;

    // Create the list of channels
    // (Each of these will be added 
    //  to every measurement)
    std::vector< HistFactory::Channel > channel_list;

    // Fill the list of channels
    for( unsigned int i = 0; i < xml_channel_files.size(); ++i ) {
      std::string channel_xml = xml_channel_files.at(i);
      cxcoutIHF << "Parsing Channel: " << channel_xml << std::endl;
      HistFactory::Channel channel =  ParseChannelXMLFile( channel_xml );

      // Get the histograms for the channel
      //collector.CollectHistograms( channel );
      //channel.CollectHistograms();
      channel_list.push_back( channel );
    }

    // Finally, add the channels to the measurements:
    for( unsigned int i = 0; i < measurement_list.size(); ++i) {

      HistFactory::Measurement& measurement = measurement_list.at(i);

      for( unsigned int j = 0; j < channel_list.size(); ++j ) {
	measurement.GetChannels().push_back( channel_list.at(j) );
      }
    }
  }
  catch(std::exception& e)
    {
      std::cout << e.what() << std::endl;
      throw hf_exc();
    }

  return measurement_list;

}
									     

HistFactory::Measurement ConfigParser::CreateMeasurementFromDriverNode( TXMLNode* node ) {


  HistFactory::Measurement measurement;

  // Set the default values:
  measurement.SetLumi( 1.0 );
  measurement.SetLumiRelErr( .10 );
  measurement.SetBinLow( 0 );
  measurement.SetBinHigh( 1 );
  measurement.SetExportOnly( false );

  cxcoutIHF << "Creating new measurement: " << std::endl;

  // First, get the attributes of the node
  TListIter attribIt = node->GetAttributes();
  TXMLAttr* curAttr = 0;
  while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

    if( curAttr->GetName() == TString( "" ) ) {
      cxcoutEHF << "Found XML attribute in Measurement with no name"  << std::endl;
      // ADD Output Here
      throw hf_exc();
    }
    else if( curAttr->GetName() == TString( "Name" ) ) {
      //rowTitle=curAttr->GetValue();
      measurement.SetName(  curAttr->GetValue() );
      //measurement.OutputFileName = outputFileNamePrefix+"_"+rowTitle+".root";
    }
    else if( curAttr->GetName() == TString( "Lumi" ) ) {
      measurement.SetLumi( atof(curAttr->GetValue()) );
    }
    else if( curAttr->GetName() == TString( "LumiRelErr" ) ) {
      measurement.SetLumiRelErr( atof(curAttr->GetValue()) );
    }
    else if( curAttr->GetName() == TString( "BinLow" ) ) {
      measurement.SetBinLow( atoi(curAttr->GetValue()) );
    }
    else if( curAttr->GetName() == TString( "BinHigh" ) ) {
      measurement.SetBinHigh( atoi(curAttr->GetValue()) );
    }
    else if( curAttr->GetName() == TString( "Mode" ) ) {
      cout <<"\n INFO: Mode attribute is deprecated, will ignore\n"<<endl;
    }
    else if( curAttr->GetName() == TString( "ExportOnly" ) ) {
      measurement.SetExportOnly( CheckTrueFalse(curAttr->GetValue(),"Measurement") );
    }

    else {
      cxcoutEHF << "Found unknown XML attribute in Measurement: " << curAttr->GetName()
		<< std::endl;
      throw hf_exc();
    }

  } // End Loop over attributes


  // Then, get the properties of the children nodes
  TXMLNode* child = node->GetChildren();
  while( child != 0 ) {
  
    if( child->GetNodeName() == TString( "" ) ) {
      cxcoutEHF << "Found XML child node of Measurement with no name"  << std::endl;
      throw hf_exc();
    }

    else if( child->GetNodeName() == TString( "POI" ) ) {
      if( child->GetText() == NULL ) {
	cxcoutEHF << "Error: node: " << child->GetName()
		  << " has no text." << std::endl;
	throw hf_exc();
      }
      //poi// measurement.SetPOI( child->GetText() );
     AddSubStrings( measurement.GetPOIList(), child->GetText() );
    }

    else if( child->GetNodeName() == TString( "ParamSetting" ) ) {
      TListIter paramIt = child->GetAttributes();
      TXMLAttr* curParam = 0;
      while( ( curParam = dynamic_cast< TXMLAttr* >( paramIt() ) ) != 0 ) {

	if( curParam->GetName() == TString( "" ) ) {
	  cxcoutEHF << "Error: Found tag attribute with no name in ParamSetting" << std::endl;
	  throw hf_exc();
	}
	else if( curParam->GetName() == TString( "Const" ) ) {
	  if(curParam->GetValue()==TString("True")){
	    // Fix here...?
	    if( child->GetText() == NULL ) {
	      cxcoutEHF << "Error: node: " << child->GetName()
			<< " has no text." << std::endl;
	      throw hf_exc();
	    }
	    AddSubStrings( measurement.GetConstantParams(), child->GetText() );
	  }
	}
	else if( curParam->GetName() == TString( "Val" ) ) {
	  double val = atof(curParam->GetValue());
	  if( child->GetText() == NULL ) {
	    cxcoutEHF << "Error: node: " << child->GetName()
		      << " has no text." << std::endl;
	    throw hf_exc();
	  }
	  std::vector<std::string> child_nodes = GetChildrenFromString(child->GetText());
	  for(unsigned int i = 0; i < child_nodes.size(); ++i) {
	    measurement.SetParamValue( child_nodes.at(i), val);
	  }
	  // AddStringValPairToMap( measurement.GetParamValues(), val, child->GetText() );
	}
	else {
	  cxcoutEHF << "Found tag attribute with unknown name in ParamSetting: "
		    << curAttr->GetName() << std::endl;
	  throw hf_exc();
	}
      }
    }

    else if( child->GetNodeName() == TString( "Asimov" ) ) {

      //std::string name;
      //std::map<string, double> fixedParams;

      // Now, create and configure an asimov object
      // and add it to the measurement
      RooStats::HistFactory::Asimov asimov;
      std::string ParamFixString;

      // Loop over attributes
      attribIt = child->GetAttributes();
      curAttr = 0;
      while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {
	
	if( curAttr->GetName() == TString( "" ) ) {
	  cxcoutEHF << "Error: Found tag attribute with no name in ConstraintTerm" << std::endl;
	  throw hf_exc();
	}

	else if( curAttr->GetName() == TString( "Name" ) ) {
	  std::string name = curAttr->GetValue();
	  asimov.SetName( name );
	}

	else if( curAttr->GetName() == TString( "FixParams" ) ) {
	  ParamFixString = curAttr->GetValue();
	  //std::map<std::string, double> fixedParams = ExtractParamMapFromString(FixParamList);
	  //asimov.GetFixedParams() = fixedParams;
	}

	else {
	  cxcoutEHF << "Found tag attribute with unknown name in ConstraintTerm: "
		    << curAttr->GetName() << std::endl;
	  throw hf_exc();
	}

      }

      // Add any parameters to the asimov dataset
      // to be fixed during the fitting and dataset generation
      if( ParamFixString=="" ) {
	cxcoutWHF << "Warning: Asimov Dataset with name: " << asimov.GetName()
		  << " added, but no parameters are set to be fixed" << std::endl;
      }
      else {
	AddParamsToAsimov( asimov, ParamFixString );
      }
      
      measurement.AddAsimovDataset( asimov );

    }

    else if( child->GetNodeName() == TString( "ConstraintTerm" ) ) {
      vector<string> syst; 
      string type = ""; 
      double rel = 0;

      map<string,double> gammaSyst;
      map<string,double> uniformSyst;
      map<string,double> logNormSyst;

      // Get the list of parameters in this tag:
      if( child->GetText() == NULL ) {
	cxcoutEHF << "Error: node: " << child->GetName()
		  << " has no text." << std::endl;
	throw hf_exc();
      }
      AddSubStrings(syst, child->GetText());

      // Now, loop over this tag's attributes
      attribIt = child->GetAttributes();
      curAttr = 0;
      while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

	if( curAttr->GetName() == TString( "" ) ) {
	  cxcoutEHF << "Error: Found tag attribute with no name in ConstraintTerm" << std::endl;
	  throw hf_exc();
	}

	else if( curAttr->GetName() == TString( "Type" ) ) {
	  type = curAttr->GetValue();
	}

	else if( curAttr->GetName() == TString( "RelativeUncertainty" ) ) {
	  rel = atof(curAttr->GetValue());
	}

	else {
	  cxcoutEHF << "Found tag attribute with unknown name in ConstraintTerm: "
		    << curAttr->GetName() << std::endl;
	  throw hf_exc();
	}

      } // End Loop over tag attributes


      // Now, fill the maps, depending on the type:

      // Check that the type is in the correct form:
      if( ! (type=="Gamma"     || type=="Uniform" || 
	     type=="LogNormal" || type=="NoConstraint") ) {
	cxcoutEHF << "Error: Encountered unknown type for ConstraintTerm: " << type << std::endl;
	throw hf_exc();
      }

      if (type=="Gamma" && rel!=0) {
	for (vector<string>::const_iterator it=syst.begin(); it!=syst.end(); ++it) {
	  // Fix Here...?
	  measurement.GetGammaSyst()[(*it).c_str()] = rel;
	}
      }
	
      if (type=="Uniform" && rel!=0) {
	for (vector<string>::const_iterator it=syst.begin(); it!=syst.end(); ++it) {
	  // Fix Here...?
	  measurement.GetUniformSyst()[(*it).c_str()] = rel;
	}
      }
	
      if (type=="LogNormal" && rel!=0) {
	for (vector<string>::const_iterator it=syst.begin(); it!=syst.end(); ++it) {
	  // Fix Here...?
	  measurement.GetLogNormSyst()[(*it).c_str()] = rel;
	}
      }
	
      if (type=="NoConstraint") {
	for (vector<string>::const_iterator it=syst.begin(); it!=syst.end(); ++it) {
	  // Fix Here...?
	  measurement.GetNoSyst()[(*it).c_str()] = 1.0; // MB : dummy value
	}
      }
    } // End adding of Constraint terms


    else if( IsAcceptableNode( child ) ) { ; }

    else {
    cxcoutEHF << "Found XML child of Measurement with unknown name: " << child->GetNodeName()
		<< std::endl;
      throw hf_exc();
    }

    child = child->GetNextNode();
  }

  measurement.PrintTree(oocoutI(static_cast<TObject*>(nullptr), HistFactory));

  return measurement;

}



HistFactory::Channel ConfigParser::ParseChannelXMLFile( string filen ) {

  /*
  TString lumiStr;
  lumiStr+=lumi;
  lumiStr.ReplaceAll(' ', TString());
  */

  cxcoutIHF << "Parsing file: " << filen ;

  TDOMParser xmlparser;

  // reading in the file and parse by DOM
  Int_t parseError = xmlparser.ParseFile( filen.c_str() );
  if( parseError ) { 
    cxcoutEHF << "Loading of xml document \"" << filen
	      << "\" failed" << std::endl;
    throw hf_exc();
  } 

  TXMLDocument* xmldoc = xmlparser.GetXMLDocument();
  TXMLNode* rootNode = xmldoc->GetRootNode();

  // Check that is is a CHANNEL based on the DOCTYPE

  if( rootNode->GetNodeName() != TString( "Channel" ) ){
    cxcoutEHF << "Error: In parsing a Channel XML, "
	      << "Encounterd XML with DOCTYPE: " << rootNode->GetNodeName() 
	      << std::endl;
    cxcoutEHF << " DOCTYPE for channels must be 'Channel' "
	      << " Check that your XML is properly written" << std::endl;
    throw hf_exc();
  }

  // Now, create the channel, 
  // configure it based on the XML
  // and return it

  HistFactory::Channel channel;

  // Set the default values:
  channel.SetInputFile( "" );
  channel.SetHistoPath( "" );

  // Walk through the root node and
  // get its attributes
  TListIter attribIt = rootNode->GetAttributes();
  TXMLAttr* curAttr = 0;
  while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

    // Get the Name, Val of this node
    TString attrName    = curAttr->GetName();
    std::string attrVal = curAttr->GetValue();

    if( attrName == TString( "" ) ) {
      cxcoutEHF << " Error: Attribute for 'Channel' with no name found" << std::endl;
      throw hf_exc();
    }

    else if( attrName == TString( "Name" ) ) {
      channel.SetName( attrVal );
      cxcoutIHF << " : creating a channel named " << channel.GetName() << std::endl;
    }

    else if( attrName == TString( "InputFile" ) ) {
      cxcoutIHF << "Setting InputFile for this channel: " << attrVal << std::endl;
      channel.SetInputFile( attrVal );
      // Set the current (cached) value
      m_currentInputFile = attrVal;        
    }

    else if( curAttr->GetName() == TString( "HistoPath" ) ) {
      cxcoutIHF << "Setting HistoPath for this channel: " << attrVal << std::endl;
      // Set the current (cached) value
      channel.SetHistoPath( attrVal );
      m_currentHistoPath = attrVal;  
    }

    else if( curAttr->GetName() == TString( "HistoName" ) ) {
      // Changed This:
      cxcoutEHF << "Use of HistoName in Channel is deprecated" << std::endl;
      cxcoutEHF << "This will be ignored" << std::endl;
    }
      
    else {
      cxcoutEHF << " Error: Unknown attribute for 'Channel' encountered: "
		<< attrName << std::endl;
      throw hf_exc();
    }

  } // End loop over the channel tag's attributes
    
  // Check that the channel was properly initiated:
  
  if( channel.GetName() == "" ) {
    cxcoutEHF << "Error: Channel created with no name" << std::endl;
    throw hf_exc();
  }

  m_currentChannel = channel.GetName();

  // Loop over the children nodes in the XML file
  // and configure the channel based on them

  TXMLNode* node = rootNode->GetChildren();

  bool firstData=true;

  while( node != 0 ) {

    // Restore the Channel-Wide Defaults
    m_currentInputFile = channel.GetInputFile();
    m_currentHistoPath = channel.GetHistoPath();
    
    if( node->GetNodeName() == TString( "" ) ) {
      cxcoutEHF << "Error: Encountered node in Channel with no name" << std::endl;
      throw hf_exc();
    }

    else if( node->GetNodeName() == TString( "Data" ) ) {
      if( firstData ) {
	RooStats::HistFactory::Data data = CreateDataElement(node);
	if( data.GetName() != "" ) {
	  cxcoutEHF << "Error: You can only rename the datasets of additional data sets.  "
		    << "  Remove the 'Name=" << data.GetName() << "' tag"
		    << " from channel: " << channel.GetName() << std::endl;
	  throw hf_exc();
	}
	channel.SetData( data );
	firstData=false;
      }
      else {
	channel.AddAdditionalData( CreateDataElement(node) );
      }
    }

    else if( node->GetNodeName() == TString( "StatErrorConfig" ) ) {
      channel.SetStatErrorConfig( CreateStatErrorConfigElement(node) );
    }

    else if( node->GetNodeName() == TString( "Sample" ) ) {
      channel.GetSamples().push_back( CreateSampleElement(node) );
    }

    else if( IsAcceptableNode( node ) ) { ; }

    else {
      cxcoutEHF << "Error: Encountered node in Channel with unknown name: " << node->GetNodeName() << std::endl;
      throw hf_exc();
    }

    node = node->GetNextNode();

  } // End loop over tags in this channel

  cxcoutIHF << "Created Channel: " << std::endl;
  channel.Print(oocoutI(static_cast<TObject*>(nullptr), HistFactory));

  return channel;

}



HistFactory::Data ConfigParser::CreateDataElement( TXMLNode* node ) {

  cxcoutIHF << "Creating Data Element" << std::endl;

    HistFactory::Data data;

    // Set the default values
    data.SetInputFile( m_currentInputFile );
    data.SetHistoPath( m_currentHistoPath );
    //data.HistoName = m_currentHistoName;

    // Now, set the attributes
    TListIter attribIt = node->GetAttributes();
    TXMLAttr* curAttr = 0;
    while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

      // Get the Name, Val of this node
      TString attrName    = curAttr->GetName();
      std::string attrVal = curAttr->GetValue();

      if( attrName == TString( "" ) ) {
	cxcoutEHF << " Error: Attribute for 'Data' with no name found" << std::endl;
	throw hf_exc();
      }

      else if( attrName == TString( "Name" ) ) {
	data.SetName( attrVal );
      }

      else if( attrName == TString( "InputFile" ) ) {
	data.SetInputFile( attrVal );
      }

      else if( attrName == TString( "HistoName" ) ) {
	data.SetHistoName( attrVal );
      }

      else if( attrName == TString( "HistoPath" ) ) {
	data.SetHistoPath( attrVal );
      }

    else if( IsAcceptableNode( node ) ) { ; }

      else {
	cxcoutEHF << " Error: Unknown attribute for 'Data' encountered: " << attrName << std::endl;
	throw hf_exc();
      }

    }

    // Check the properties of the data node:
    if( data.GetInputFile() == "" ) {
      cxcoutEHF << "Error: Data Node has no InputFile" << std::endl;
      throw hf_exc();
    }
    if( data.GetHistoName() == "" ) {
      cxcoutEHF << "Error: Data Node has no HistoName" << std::endl;
      throw hf_exc();
    }

    cxcoutIHF << "Created Data Node with"
	      << " InputFile: " << data.GetInputFile()
	      << " HistoName: " << data.GetHistoName()
	      << " HistoPath: " << data.GetHistoPath()
        << (data.GetName() != "" ? " Name: " : "") << data.GetName() << std::endl;

    // data.hist = GetHisto(data.FileName, data.HistoPath, data.HistoName);

    return data;
}



HistFactory::StatErrorConfig ConfigParser::CreateStatErrorConfigElement( TXMLNode* node ) {

  cxcoutIHF << "Creating StatErrorConfig Element" << std::endl;
  
  HistFactory::StatErrorConfig config;

  // Setup default values:
  config.SetConstraintType( Constraint::Poisson );
  config.SetRelErrorThreshold( 0.05 ); // 5%

  // Loop over the node's attributes
  TListIter attribIt = node->GetAttributes();
  TXMLAttr* curAttr = 0;
  while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

    // Get the Name, Val of this node
    TString attrName    = curAttr->GetName();
    std::string attrVal = curAttr->GetValue();
    
    if( attrName == TString( "RelErrorThreshold" ) ) {
      config.SetRelErrorThreshold( atof(attrVal.c_str()) );
    }
    
    if( attrName == TString( "ConstraintType" ) ) {
      // Allowable Values:  Gaussian

      if( attrVal == "" ) {
	cxcoutEHF << "Error: Bad Value for StatErrorConfig Constraint Type Found" << std::endl;
	throw hf_exc();
      }

      else if( attrVal=="Gaussian" || attrVal=="Gauss"  ) {  
	config.SetConstraintType( Constraint::Gaussian );
      }

      else if( attrVal=="Poisson" || attrVal=="Pois"  ) {
	config.SetConstraintType( Constraint::Poisson );
      }

      else if( IsAcceptableNode( node ) ) { ; }

      else {
	cout << "Invalid Stat Constraint Type: " << curAttr->GetValue() << endl;
	throw hf_exc();
      }
    }
  } // End: Loop Over Attributes

  cxcoutIHF << "Created StatErrorConfig Element with"
	    << " Constraint type: " << config.GetConstraintType()
	    << " RelError Threshold: " << config.GetRelErrorThreshold()
	    << std::endl;

  return config;

}


HistFactory::Sample ConfigParser::CreateSampleElement( TXMLNode* node ) {

  cxcoutIHF << "Creating Sample Element" << std::endl;

  HistFactory::Sample sample;

  // Set the default values
  sample.SetInputFile( m_currentInputFile );
  sample.SetHistoPath( m_currentHistoPath );
  sample.SetChannelName( m_currentChannel );
  sample.SetNormalizeByTheory( true );
  //sample.HistoName = m_currentHistoName;

  // Now, set the attributes

  TListIter attribIt = node->GetAttributes();
  TXMLAttr* curAttr = 0;
  while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

    // Get the Name, Val of this node
    TString attrName    = curAttr->GetName();
    std::string attrVal = curAttr->GetValue();

    if( attrName == TString( "" ) ) {
      cxcoutEHF << " Error: Attribute for 'Sample' with no name found" << std::endl;
      throw hf_exc();
    }

    else if( attrName == TString( "Name" ) ) {
      sample.SetName( attrVal );
    }

    else if( attrName == TString( "InputFile" ) ) {
      sample.SetInputFile( attrVal );
      m_currentInputFile = attrVal;
    }

    else if( attrName == TString( "HistoName" ) ) {
      sample.SetHistoName( attrVal );
    }

    else if( attrName == TString( "HistoPath" ) ) {
      sample.SetHistoPath( attrVal );
      m_currentHistoPath = attrVal;
    }

    else if( attrName == TString( "NormalizeByTheory" ) ) {
      sample.SetNormalizeByTheory( CheckTrueFalse(attrVal,"Sample") );
      /*
	if( attrVal == "" ) {
	cxcoutEHF << "Error: Attribute 'NormalizeByTheory' in Sample has no value" << std::endl;
	throw hf_exc();
	}
	else if ( attrVal == "True"  || attrVal == "true"  )   sample.NormalizeByTheory = true;
	else if ( attrVal == "False" || attrVal == "false" )   sample.NormalizeByTheory = false;
	else {
	cxcoutEHF << "Error: Attribute 'NormalizeByTheory' in Sample has unknown value: " << attrVal <<  std::endl;
	std::cout << "Value must be 'True' or 'False' " <<  std::endl;
	throw hf_exc();
	}
      */
    }

    else {
      cxcoutEHF << " Error: Unknown attribute for 'Sample' encountered: " << attrName << std::endl;
      throw hf_exc();
    }
  }

  // Quickly check the properties of the Sample Node
  if( sample.GetName() == "" ) {
    cxcoutEHF << "Error: Sample Node has no Name" << std::endl;
    throw hf_exc();
  }
  if( sample.GetInputFile() == "" ) {
    cxcoutEHF << "Error: Sample Node has no InputFile" << std::endl;
    throw hf_exc();
  }
  if( sample.GetHistoName() == "" ) {
    cxcoutEHF << "Error: Sample Node has no HistoName" << std::endl;
    throw hf_exc();
  }


  // Now, loop over the children and add the systematics

  TXMLNode* child = node->GetChildren();

  while( child != 0 ) {
      
    if( child->GetNodeName() == TString( "" ) ) {
      cxcoutEHF << "Error: Encountered node in Sample with no name" << std::endl;
      throw hf_exc();
    }

    else if( child->GetNodeName() == TString( "NormFactor" ) ) {
      sample.GetNormFactorList().push_back( MakeNormFactor( child ) );
    }

    else if( child->GetNodeName() == TString( "OverallSys" ) ) {
      sample.GetOverallSysList().push_back( MakeOverallSys( child ) );
    }

    else if( child->GetNodeName() == TString( "HistoSys" ) ) {
      sample.GetHistoSysList().push_back( MakeHistoSys( child ) );
    }

    else if( child->GetNodeName() == TString( "HistoFactor" ) ) {
      cxcoutEHF << "WARNING: HistoFactor not yet supported" << std::endl;
      //sample.GetHistoFactorList().push_back( MakeHistoFactor( child ) );
    }

    else if( child->GetNodeName() == TString( "ShapeSys" ) ) {
      sample.GetShapeSysList().push_back( MakeShapeSys( child ) );
    }

    else if( child->GetNodeName() == TString( "ShapeFactor" ) ) {
      sample.GetShapeFactorList().push_back( MakeShapeFactor( child ) );
    }

    else if( child->GetNodeName() == TString( "StatError" ) ) {
      sample.SetStatError( ActivateStatError(child) );
    }

    else if( IsAcceptableNode( child ) ) { ; }

    else {
      cxcoutEHF << "Error: Encountered node in Sample with unknown name: " << child->GetNodeName() << std::endl;
      throw hf_exc();
    }

    child=child->GetNextNode();
  }

  cxcoutIHF << "Created Sample Node with"
	    << " Name: " << sample.GetName()
	    << " InputFile: " << sample.GetInputFile()
	    << " HistoName: " << sample.GetHistoName()
	    << " HistoPath: " << sample.GetHistoPath()
	    << std::endl;

  // sample.hist = GetHisto(sample.FileName, sample.HistoPath, sample.HistoName);

  return sample;
}


HistFactory::NormFactor ConfigParser::MakeNormFactor( TXMLNode* node ) {

  cxcoutIHF << "Making NormFactor:" << std::endl;
  
  HistFactory::NormFactor norm;

  TListIter attribIt = node->GetAttributes();
  TXMLAttr* curAttr = 0;
  while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

    // Get the Name, Val of this node
    TString attrName    = curAttr->GetName();
    std::string attrVal = curAttr->GetValue();

    if( attrName == TString( "" ) ){
      cxcoutEHF << "Error: Encountered Element in NormFactor with no name" << std::endl;
      throw hf_exc();
    }

    else if( curAttr->GetName() == TString( "Name" ) ) {
      norm.SetName( attrVal );
    }
    else if( curAttr->GetName() == TString( "Val" ) ) {
      norm.SetVal( atof(attrVal.c_str()) );
    }
    else if( curAttr->GetName() == TString( "Low" ) ) {
      norm.SetLow( atof(attrVal.c_str()) );
    }
    else if( curAttr->GetName() == TString( "High" ) ) {
      norm.SetHigh( atof(attrVal.c_str()) );
    }
    else if( curAttr->GetName() == TString( "Const" ) ) {
      norm.SetConst( CheckTrueFalse(attrVal,"NormFactor") );
    }

    else {
      cxcoutEHF << "Error: Encountered Element in NormFactor with unknown name: "
		<< attrName << std::endl;
      throw hf_exc();
    }

  } // End loop over properties

  if( norm.GetName() == "" ) {
    cxcoutEHF << "Error: NormFactor Node has no Name" << std::endl;
    throw hf_exc();
  }

  if( norm.GetLow() >= norm.GetHigh() ) {
    cxcoutEHF << "Error: NormFactor: " << norm.GetName()
	      << " has lower limit >= its upper limit: " 
	      << " Lower: " << norm.GetLow()
	      << " Upper: " << norm.GetHigh()
	      << ". Please Fix" << std::endl;
    throw hf_exc();
  }
  if( norm.GetVal() > norm.GetHigh() || norm.GetVal() < norm.GetLow() ) {
    cxcoutEHF << "Error: NormFactor: " << norm.GetName()
	      << " has initial value not within its range: "
	      << " Val: " << norm.GetVal()
	      << " Lower: " << norm.GetLow()
	      << " Upper: " << norm.GetHigh()
	      << ". Please Fix" << std::endl;
    throw hf_exc();
  }

  norm.Print(oocoutI(static_cast<TObject*>(nullptr), HistFactory));

  return norm;

}

HistFactory::HistoFactor ConfigParser::MakeHistoFactor( TXMLNode* /*node*/ ) {

  cxcoutIHF << "Making HistoFactor" << std::endl;

  HistFactory::HistoFactor dummy;

  dummy.SetInputFileLow( m_currentInputFile );
  dummy.SetHistoPathLow( m_currentHistoPath );

  dummy.SetInputFileHigh( m_currentInputFile );
  dummy.SetHistoPathHigh( m_currentHistoPath );

  cxcoutIHF << "Made HistoFactor" << std::endl;

  return dummy;

}


HistFactory::HistoSys ConfigParser::MakeHistoSys( TXMLNode* node ) {

  cxcoutIHF << "Making HistoSys:" << std::endl;

  HistFactory::HistoSys histoSys;

  // Set Default values
  histoSys.SetInputFileLow( m_currentInputFile );
  histoSys.SetHistoPathLow( m_currentHistoPath );

  histoSys.SetInputFileHigh( m_currentInputFile );
  histoSys.SetHistoPathHigh( m_currentHistoPath );

  TListIter attribIt = node->GetAttributes();
  TXMLAttr* curAttr = 0;
  /*
  string Name, histoPathHigh, histoPathLow, 
    histoNameLow, histoNameHigh, inputFileHigh, inputFileLow;
  inputFileLow=inputFileName; inputFileHigh=inputFileName;
  histoPathLow=histoPathName; histoPathHigh=histoPathName;
  histoNameLow=histoName; histoNameHigh=histoName;
  */

  while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

    // Get the Name, Val of this node
    TString attrName    = curAttr->GetName();
    std::string attrVal = curAttr->GetValue();
	
    if( attrName == TString( "" ) ){
      cxcoutEHF << "Error: Encountered Element in HistoSys with no name" << std::endl;
      throw hf_exc();
    }

    else if( curAttr->GetName() == TString( "Name" ) ) {
      histoSys.SetName( attrVal );
    }

    else if( curAttr->GetName() == TString( "HistoFileHigh" ) ) {
      histoSys.SetInputFileHigh( attrVal );
    }
    else if( curAttr->GetName() == TString( "HistoPathHigh" ) ) {
      histoSys.SetHistoPathHigh( attrVal );
    }
    else if( curAttr->GetName() == TString( "HistoNameHigh" ) ) {
      histoSys.SetHistoNameHigh( attrVal );
    }

    else if( curAttr->GetName() == TString( "HistoFileLow" ) ) {
      histoSys.SetInputFileLow( attrVal );
    }
    else if( curAttr->GetName() == TString( "HistoPathLow" ) ) {
      histoSys.SetHistoPathLow( attrVal );
    }
    else if( curAttr->GetName() == TString( "HistoNameLow" ) ) {
      histoSys.SetHistoNameLow( attrVal );
    }

    else {
      cxcoutEHF << "Error: Encountered Element in HistoSys with unknown name: "
		<< attrName << std::endl;
      throw hf_exc();
    }

  } // End loop over properties


  if( histoSys.GetName() == "" ) {
    cxcoutEHF << "Error: HistoSys Node has no Name" << std::endl;
    throw hf_exc();
  }
  if( histoSys.GetInputFileHigh() == "" ) {
    cxcoutEHF << "Error: HistoSysSample Node has no InputFileHigh" << std::endl;
    throw hf_exc();
  }
  if( histoSys.GetInputFileLow() == "" ) {
    cxcoutEHF << "Error: HistoSysSample Node has no InputFileLow" << std::endl;
    throw hf_exc();
  }
  if( histoSys.GetHistoNameHigh() == "" ) {
    cxcoutEHF << "Error: HistoSysSample Node has no HistoNameHigh" << std::endl;
    throw hf_exc();
  }
  if( histoSys.GetHistoNameLow() == "" ) {
    cxcoutEHF << "Error: HistoSysSample Node has no HistoNameLow" << std::endl;
    throw hf_exc();
  }


  histoSys.Print(oocoutI(static_cast<TObject*>(nullptr), HistFactory));

  return histoSys;

}


HistFactory::OverallSys ConfigParser::MakeOverallSys( TXMLNode* node ) {

  cxcoutIHF << "Making OverallSys:" << std::endl;

  HistFactory::OverallSys overallSys;

  TListIter attribIt = node->GetAttributes();
  TXMLAttr* curAttr = 0;
  while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

    // Get the Name, Val of this node
    TString attrName    = curAttr->GetName();
    std::string attrVal = curAttr->GetValue();

    if( attrName == TString( "" ) ){
      cxcoutEHF << "Error: Encountered Element in OverallSys with no name" << std::endl;
      throw hf_exc();
    }

    else if( attrName == TString( "Name" ) ) {
      overallSys.SetName( attrVal );
    }
    else if( attrName == TString( "High" ) ) {
      overallSys.SetHigh( atof(attrVal.c_str()) );
    }
    else if( attrName == TString( "Low" ) ) {
      overallSys.SetLow( atof(attrVal.c_str()) );
    }

    else {
      cxcoutEHF << "Error: Encountered Element in OverallSys with unknown name: "
		<< attrName << std::endl;
      throw hf_exc();
    }

  }

  if( overallSys.GetName() == "" ) {
    cxcoutEHF << "Error: Encountered OverallSys with no name" << std::endl;
    throw hf_exc();
  }
  

  overallSys.Print(oocoutI(static_cast<TObject*>(nullptr), HistFactory));

  return overallSys;

}


HistFactory::ShapeFactor ConfigParser::MakeShapeFactor( TXMLNode* node ) {

  cxcoutIHF << "Making ShapeFactor" << std::endl;

  HistFactory::ShapeFactor shapeFactor;

  TListIter attribIt = node->GetAttributes();
  TXMLAttr* curAttr = 0;

  // A Shape Factor may or may not include an initial shape
  // This will be set by strings pointing to a histogram
  // If we don't see a 'HistoName' attribute, we assume
  // that an initial shape is not being set
  std::string ShapeInputFile = m_currentInputFile;
  std::string ShapeInputPath = m_currentHistoPath;

  while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

    // Get the Name, Val of this node
    TString attrName    = curAttr->GetName();
    std::string attrVal = curAttr->GetValue();

    if( attrName == TString( "" ) ){
      cxcoutEHF << "Error: Encountered Element in ShapeFactor with no name" << std::endl;
      throw hf_exc();
    }

    else if( attrName == TString( "Name" ) ) {
      shapeFactor.SetName( attrVal );
    }
    else if( attrName == TString( "Const" ) ) {
      shapeFactor.SetConstant( CheckTrueFalse(attrVal, "ShapeFactor" ) );
    }
    
    else if( attrName == TString( "HistoName" ) ) {
      shapeFactor.SetHistoName( attrVal );
    }
    
    else if( attrName == TString( "InputFile" ) ) {
      ShapeInputFile = attrVal;
    }
    
    else if( attrName == TString( "HistoPath" ) ) {
      ShapeInputPath = attrVal;
    }

    else {
      cxcoutEHF << "Error: Encountered Element in ShapeFactor with unknown name: "
		<< attrName << std::endl;
      throw hf_exc();
    }

  }

  if( shapeFactor.GetName() == "" ) {
    cxcoutEHF << "Error: Encountered ShapeFactor with no name" << std::endl;
    throw hf_exc();
  }

  // Set the Histogram name, path, and file
  // if an InitialHist is set
  if( shapeFactor.HasInitialShape() ) {
    if( shapeFactor.GetHistoName() == "" ) {
      cxcoutEHF << "Error: ShapeFactor: " << shapeFactor.GetName()
		<< " is configured to have an initial shape, but "
		<< "its histogram doesn't have a name"
		<< std::endl;
      throw hf_exc();
    }
    shapeFactor.SetHistoPath( ShapeInputPath );
    shapeFactor.SetInputFile( ShapeInputFile );
  }
  
  shapeFactor.Print(oocoutI(static_cast<TObject*>(nullptr), HistFactory));

  return shapeFactor;

}


HistFactory::ShapeSys ConfigParser::MakeShapeSys( TXMLNode* node ) {

  cxcoutIHF << "Making ShapeSys" << std::endl;

  HistFactory::ShapeSys shapeSys;

  // Set the default values
  shapeSys.SetConstraintType( Constraint::Gaussian );
  shapeSys.SetInputFile( m_currentInputFile );
  shapeSys.SetHistoPath( m_currentHistoPath );


  TListIter attribIt = node->GetAttributes();
  TXMLAttr* curAttr = 0;
  //EstimateSummary::ConstraintType ConstraintType = EstimateSummary::Gaussian; //"Gaussian";

  while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {


    // Get the Name, Val of this node
    TString attrName    = curAttr->GetName();
    std::string attrVal = curAttr->GetValue();

    if( attrName == TString( "" ) ){
      cxcoutEHF << "Error: Encountered Element in ShapeSys with no name" << std::endl;
      throw hf_exc();
    }

    else if( attrName == TString( "Name" ) ) {
      shapeSys.SetName( attrVal );
    }

    else if( attrName == TString( "HistoName" ) ) {
      shapeSys.SetHistoName( attrVal );
    }

    else if( attrName == TString( "HistoPath" ) ) {
      shapeSys.SetHistoPath( attrVal );
    }

    else if( attrName == TString( "InputFile" ) ) {
      shapeSys.SetInputFile( attrVal );
    }

    else if( attrName == TString( "ConstraintType" ) ) {
      if( attrVal=="" ) {
	cxcoutEHF << "Error: ShapeSys Constraint type is empty" << std::endl;
	throw hf_exc();
      }
      else if( attrVal=="Gaussian" || attrVal=="Gauss" ) {
	shapeSys.SetConstraintType( Constraint::Gaussian );
      }
      else if( attrVal=="Poisson"  || attrVal=="Pois"  ) {
	shapeSys.SetConstraintType( Constraint::Poisson );
      }
      else {
	cout << "Error: Encountered unknown ShapeSys Constraint type: " << attrVal << endl;
	throw hf_exc();
      }
    }

    else {
      cxcoutEHF << "Error: Encountered Element in ShapeSys with unknown name: "
		<< attrName << std::endl;
      throw hf_exc();
    }

  } // End loop over attributes


  if( shapeSys.GetName() == "" ) {
    cxcoutEHF << "Error: Encountered ShapeSys with no Name" << std::endl;
    throw hf_exc();
  }
  if( shapeSys.GetInputFile() == "" ) {
    cxcoutEHF << "Error: Encountered ShapeSys with no InputFile" << std::endl;
    throw hf_exc();
  }
  if( shapeSys.GetHistoName() == "" ) {
    cxcoutEHF << "Error: Encountered ShapeSys with no HistoName" << std::endl;
    throw hf_exc();
  }

  shapeSys.Print(oocoutI(static_cast<TObject*>(nullptr), HistFactory));

  return shapeSys;

}


HistFactory::StatError ConfigParser::ActivateStatError( TXMLNode* node ) {
	
  cxcoutIHF << "Activating StatError" << std::endl;

  // Set default values
  HistFactory::StatError statError;
  statError.Activate( false );
  statError.SetUseHisto( false );
  statError.SetHistoName( "" );

  // Loop over the node's attributes
  TListIter attribIt = node->GetAttributes();
  TXMLAttr* curAttr = 0;
  while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {

    // Get the Name, Val of this node
    TString attrName    = curAttr->GetName();
    std::string attrVal = curAttr->GetValue();

    if( attrName == TString( "" ) ){
      cxcoutEHF << "Error: Encountered Element in ActivateStatError with no name" << std::endl;
      throw hf_exc();
    }
    
    else if( attrName == TString( "Activate" ) ) {
      statError.Activate( CheckTrueFalse(attrVal,"ActivateStatError") );
    }

    else if( attrName == TString( "HistoName" ) ) {
      statError.SetHistoName( attrVal );
    }

    else if( attrName == TString( "HistoPath" ) ) {
      statError.SetHistoPath( attrVal );
    }

    else if( attrName == TString( "InputFile" ) ) {
      statError.SetInputFile( attrVal );
    }
    
    else {
      cxcoutEHF << "Error: Encountered Element in ActivateStatError with unknown name: "
		<< attrName << std::endl;
      throw hf_exc();
    }

  } // End: Loop Over Attributes

  // Based on the input, determine
  // if we should use a histogram or not:
  // Logic: One turns on using a histogram
  // by setting the attribute "HistoName"
  // If this is set AND the InputFile or
  // HistoPath aren't set, we set those
  // to the current default values
  if( statError.GetHistoName() != "" ) {
    statError.SetUseHisto( true );

    // Check that a file has been set
    // (Possibly using the default)
    if( statError.GetInputFile() == "" ) {
      statError.SetInputFile( m_currentInputFile );
    }
    if( statError.GetHistoPath() == "" ) {
      statError.SetHistoPath( m_currentHistoPath );
    }

  }

  /*
  if( statError.Activate ) {
    if( statError.UseHisto ) {
    }
    else {
      statError.InputFile = "";
      statError.HistoName = "";
      statError.HistoPath = "";
    }
    }
  */

  statError.Print(oocoutI(static_cast<TObject*>(nullptr), HistFactory));

  return statError;

}


RooStats::HistFactory::PreprocessFunction ConfigParser::ParseFunctionConfig( TXMLNode* functionNode ){

  cxcoutIHF << "Parsing FunctionConfig" << std::endl;

  //std::string name, expression, dependents;
  TListIter attribIt = functionNode->GetAttributes();
  TXMLAttr* curAttr = 0;

  std::string Name = "";
  std::string Expression = "";
  std::string Dependents = "";

  // Add protection to ensure that all parts are there
  while( ( curAttr = dynamic_cast< TXMLAttr* >( attribIt() ) ) != 0 ) {
    if( curAttr->GetName() == TString( "Name" ) ) {
      Name = curAttr->GetValue();
      //func.SetName( curAttr->GetValue() );
      //name = curAttr->GetValue() ;
    }
    if( curAttr->GetName() == TString( "Expression" ) ) {
      Expression = curAttr->GetValue();
      //func.SetExpression( curAttr->GetValue() );
    }
    if( curAttr->GetName() == TString( "Dependents" ) ) {
      Dependents = curAttr->GetValue();
      //func.SetDependents( curAttr->GetValue() );
    }    
  }
  
  if( Name=="" ){
    cxcoutEHF << "Error processing PreprocessFunction: Name attribute is empty" << std::endl;
    throw hf_exc();
  }
  if( Expression=="" ){
    cxcoutEHF << "Error processing PreprocessFunction: Expression attribute is empty" << std::endl;
    throw hf_exc();
  }
  if( Dependents=="" ){
    cxcoutEHF << "Error processing PreprocessFunction: Dependents attribute is empty" << std::endl;
    throw hf_exc();
  }

  RooStats::HistFactory::PreprocessFunction func(Name, Expression, Dependents);
  
  cxcoutIHF << "Created Preprocess Function: " << func.GetCommand() << std::endl;

  //std::string command = "expr::"+func.GetName()+"('"+func.GetExpression()+"',{"+func.GetDependents()+"})";
  //func.SetCommand( command );
  // //  cout << "will pre-process this line " << ret <<endl;
  return func;

}


bool ConfigParser::IsAcceptableNode( TXMLNode* node ) {

  if( node->GetNodeName() == TString( "text" ) ) {    
    return true;
  }

  if( node->GetNodeName() == TString( "comment" ) ) {
    return true;
  }

  return false;

}


bool ConfigParser::CheckTrueFalse( std::string attrVal, std::string NodeTitle ) {

  if( attrVal == "" ) {
    cxcoutEHF << "Error: In " << NodeTitle
	      << " Expected either 'True' or 'False' but found empty" << std::endl;
    throw hf_exc();
  }
  else if ( attrVal == "True"  || attrVal == "true"  )   return true;
  else if ( attrVal == "False" || attrVal == "false" )   return false;
  else {
    cxcoutEHF << "Error: In " << NodeTitle
	      << " Expected either 'True' or 'False' but found: " << attrVal <<  std::endl;
    throw hf_exc();
  }

  return false;
  
}

//ConfigParser
