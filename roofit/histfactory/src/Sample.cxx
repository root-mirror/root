// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, George Lewis 
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////////
/** \class RooStats::HistFactory::Sample
 * \ingroup HistFactory
 */

#include "TH1.h"
#include "RooStats/HistFactory/Sample.h"
#include "RooStats/HistFactory/HistFactoryException.h"

//#include "TClass.h"

RooStats::HistFactory::Sample::Sample() : 
  fNormalizeByTheory(false), fStatErrorActivate(false), fhNominal(), fhCountingHist(0) { ; }

// copy constructor (important for python)
RooStats::HistFactory::Sample::Sample(const Sample& other) :
  fName(other.fName), fInputFile(other.fInputFile),
  fHistoName(other.fHistoName), fHistoPath(other.fHistoPath),
  fChannelName(other.fChannelName),

  fOverallSysList(other.fOverallSysList),
  fNormFactorList(other.fNormFactorList),
  fHistoSysList(other.fHistoSysList),
  fHistoFactorList(other.fHistoFactorList),
  fShapeSysList(other.fShapeSysList),
  fShapeFactorList(other.fShapeFactorList),

  fStatError(other.fStatError),
  fNormalizeByTheory(other.fNormalizeByTheory),
  fStatErrorActivate(other.fStatErrorActivate),
  fhNominal(other.fhNominal),
  fhCountingHist(0)
  { 
    if( other.fhCountingHist ) {
      SetValue( other.fhCountingHist->GetBinContent(1) );
    }else{
      fhCountingHist = NULL;
    }
  }

RooStats::HistFactory::Sample& RooStats::HistFactory::Sample::operator=(const Sample& other)
{
  fName = other.fName; fInputFile = other.fInputFile;
  fHistoName = other.fHistoName; fHistoPath = other.fHistoPath;
  fChannelName = other.fChannelName;

  fOverallSysList = other.fOverallSysList;
  fNormFactorList = other.fNormFactorList;
  fHistoSysList = other.fHistoSysList;
  fHistoFactorList = other.fHistoFactorList;
  fShapeSysList = other.fShapeSysList;
  fShapeFactorList = other.fShapeFactorList;

  fStatError = other.fStatError;
  fNormalizeByTheory = other.fNormalizeByTheory;
  fStatErrorActivate = other.fStatErrorActivate;
  fhNominal = other.fhNominal;

  if (fhCountingHist)
    delete fhCountingHist;

  if( other.fhCountingHist ) {
    SetValue( other.fhCountingHist->GetBinContent(1) );
  } else {
    fhCountingHist = NULL;
  }

  return *this;
}


RooStats::HistFactory::Sample::Sample(std::string SampName, std::string SampHistoName, std::string SampInputFile, std::string SampHistoPath) : 
  fName( SampName ),   fInputFile( SampInputFile), 
  fHistoName( SampHistoName ), fHistoPath( SampHistoPath ),
  fNormalizeByTheory(true), fStatErrorActivate(false), fhNominal(),
  fhCountingHist(0) { ; }

RooStats::HistFactory::Sample::Sample(std::string SampName) : 
  fName( SampName ),   fInputFile( "" ), 
  fHistoName( "" ), fHistoPath( "" ),
  fNormalizeByTheory(true), fStatErrorActivate(false),fhNominal(),
  fhCountingHist(0) { ; }

RooStats::HistFactory::Sample::~Sample() {
  if(fhCountingHist)
    delete fhCountingHist;
}

const TH1* RooStats::HistFactory::Sample::GetHisto() const {
  TH1* histo = (TH1*) fhNominal.GetObject();
  return histo;
}


void RooStats::HistFactory::Sample::writeToFile( std::string OutputFileName, std::string DirName ) {

  const TH1* histNominal = GetHisto();
  histNominal->Write();
  
  // Set the location of the data
  // in the output measurement
  
  fInputFile = OutputFileName;
  fHistoName = histNominal->GetName();
  fHistoPath = DirName;

  // Write this sample's StatError
  GetStatError().writeToFile( OutputFileName, DirName );

  // Must write all systematics that contain internal histograms
  // (This is not all systematics)
  for( unsigned int i = 0; i < GetHistoSysList().size(); ++i ) {
    GetHistoSysList().at(i).writeToFile( OutputFileName, DirName );
  }
  for( unsigned int i = 0; i < GetHistoFactorList().size(); ++i ) {
    GetHistoFactorList().at(i).writeToFile( OutputFileName, DirName );
  }
  for( unsigned int i = 0; i < GetShapeSysList().size(); ++i ) {
    GetShapeSysList().at(i).writeToFile( OutputFileName, DirName );
  }
  for( unsigned int i = 0; i < GetShapeFactorList().size(); ++i ) {
    GetShapeFactorList().at(i).writeToFile( OutputFileName, DirName );
  }

  return;

}


void RooStats::HistFactory::Sample::SetValue( Double_t val ) {

  // For use in a number counting measurement
  // Create a 1-bin histogram, 
  // fill it with this input value,
  // and set this Sample's histogram to that hist
  
  std::string SampleHistName = fName + "_hist";
  
  // Histogram has 1-bin (hard-coded)
  if(fhCountingHist)
    delete fhCountingHist;
  
  fhCountingHist = new TH1F( SampleHistName.c_str(), SampleHistName.c_str(), 1, 0, 1 );
  fhCountingHist->SetBinContent( 1, val );

  // Set the histogram of the internally held data
  // node of this channel to this newly created histogram
  SetHisto( fhCountingHist );

}



void RooStats::HistFactory::Sample::Print( std::ostream& stream ) const {


  stream << "\t \t Name: " << fName
	 << "\t \t Channel: " << fChannelName
	 << "\t NormalizeByTheory: " << (fNormalizeByTheory ? "True" : "False")
	 << "\t StatErrorActivate: " << (fStatErrorActivate ? "True" : "False")
	 << std::endl;  

  stream << "\t \t \t \t " 
	 << "\t InputFile: " << fInputFile
	 << "\t HistName: " << fHistoName
	 << "\t HistoPath: " << fHistoPath
	 << "\t HistoAddress: " << GetHisto()
    // << "\t Type: " << GetHisto()->ClassName()
	 << std::endl;  

  if( fStatError.GetActivate() ) {
    stream << "\t \t \t StatError Activate: " << fStatError.GetActivate()
	   << "\t InputFile: " << fInputFile
	   << "\t HistName: " << fStatError.GetHistoName()
	   << "\t HistoPath: " << fStatError.GetHistoPath()
	   << "\t HistoAddress: " << fStatError.GetErrorHist()
         << std::endl;  
  }


  /*
  stream<< " NormalizeByTheory: ";
  if(NormalizeByTheory)  stream << "True";
  else                   stream << "False";

  stream<< " StatErrorActivate: ";
  if(StatErrorActivate)  stream << "True";
  else                   stream << "False";
  */


}  

void RooStats::HistFactory::Sample::PrintXML( std::ofstream& xml ) {
  
  
  // Create the sample tag
  xml << "    <Sample Name=\"" << fName << "\" "
      << " HistoPath=\"" << fHistoPath << "\" "
      << " HistoName=\"" << fHistoName << "\" "
      << " InputFile=\"" << fInputFile << "\" "
      << " NormalizeByTheory=\"" << (fNormalizeByTheory ? std::string("True") : std::string("False"))  << "\" "
      << ">" << std::endl;


  // Print Stat Error (if necessary)
  fStatError.PrintXML( xml );
  /*
  if( fStatError.GetActivate() ) {
    xml << "      <StatError Activate=\"" << (fStatError.GetActivate() ? std::string("True") : std::string("False"))  << "\" "
	<< " InputFile=\"" << fStatError.GetInputFile() << "\" "
	<< " HistoName=\"" << fStatError.GetHistoName() << "\" "
	<< " HistoPath=\"" << fStatError.GetHistoPath() << "\" "
	<< " /> " << std::endl;
  }
  */


  // Now, print the systematics:
  for( unsigned int i = 0; i < fOverallSysList.size(); ++i ) {
    RooStats::HistFactory::OverallSys sys = fOverallSysList.at(i);
    sys.PrintXML(xml);
    /*
    xml << "      <OverallSys Name=\"" << sys.GetName() << "\" "
	<< " High=\"" << sys.GetHigh() << "\" "
	<< " Low=\""  << sys.GetLow()  << "\" "
	<< "  /> " << std::endl;
    */
  }
  for( unsigned int i = 0; i < fNormFactorList.size(); ++i ) {
    RooStats::HistFactory::NormFactor sys = fNormFactorList.at(i);
    sys.PrintXML(xml);
    /*
    xml << "      <NormFactor Name=\"" << sys.GetName() << "\" "
	<< " Val=\""   << sys.GetVal()   << "\" "
	<< " High=\""  << sys.GetHigh()  << "\" "
	<< " Low=\""   << sys.GetLow()   << "\" "
	<< " Const=\"" << (sys.GetConst() ? std::string("True") : std::string("False")) << "\" "
	<< "  /> " << std::endl;
    */
  }
  for( unsigned int i = 0; i < fHistoSysList.size(); ++i ) {
    RooStats::HistFactory::HistoSys sys = fHistoSysList.at(i);
    sys.PrintXML(xml);
    /*
    xml << "      <HistoSys Name=\"" << sys.GetName() << "\" "

	<< " InputFileLow=\""  << sys.GetInputFileLow()  << "\" "
	<< " HistoNameLow=\""  << sys.GetHistoNameLow()  << "\" "
	<< " HistoPathLow=\""  << sys.GetHistoPathLow()  << "\" "

	<< " InputFileHigh=\""  << sys.GetInputFileHigh()  << "\" "
	<< " HistoNameHigh=\""  << sys.GetHistoNameHigh()  << "\" "
	<< " HistoPathHigh=\""  << sys.GetHistoPathHigh()  << "\" "
	<< "  /> " << std::endl;
    */
  }
  for( unsigned int i = 0; i < fHistoFactorList.size(); ++i ) {
    RooStats::HistFactory::HistoFactor sys = fHistoFactorList.at(i);
    sys.PrintXML(xml);
    /*
    xml << "      <HistoFactor Name=\"" << sys.GetName() << "\" "

	<< " InputFileLow=\""  << sys.GetInputFileLow()  << "\" "
	<< " HistoNameLow=\""  << sys.GetHistoNameLow()  << "\" "
	<< " HistoPathLow=\""  << sys.GetHistoPathLow()  << "\" "

	<< " InputFileHigh=\""  << sys.GetInputFileHigh()  << "\" "
	<< " HistoNameHigh=\""  << sys.GetHistoNameHigh()  << "\" "
	<< " HistoPathHigh=\""  << sys.GetHistoPathHigh()  << "\" "
	<< "  /> " << std::endl;
    */
  }
  for( unsigned int i = 0; i < fShapeSysList.size(); ++i ) {
    RooStats::HistFactory::ShapeSys sys = fShapeSysList.at(i);
    sys.PrintXML(xml);
    /*
    xml << "      <ShapeSys Name=\"" << sys.GetName() << "\" "

	<< " InputFile=\""  << sys.GetInputFile()  << "\" "
	<< " HistoName=\""  << sys.GetHistoName()  << "\" "
	<< " HistoPath=\""  << sys.GetHistoPath()  << "\" "
	<< " ConstraintType=\"" << std::string(Constraint::Name(sys.GetConstraintType())) << "\" "
	<< "  /> " << std::endl;
    */
  }
  for( unsigned int i = 0; i < fShapeFactorList.size(); ++i ) {
    RooStats::HistFactory::ShapeFactor sys = fShapeFactorList.at(i);
    sys.PrintXML(xml);
    /*
    xml << "      <ShapeFactor Name=\"" << sys.GetName() << "\" "
	<< "  /> " << std::endl;
    */
  }

  // Finally, close the tag
  xml << "    </Sample>" << std::endl;

}


#ifdef INCLUDE_RYML
#include <ryml.hpp>
#include <c4/yml/std/map.hpp>
#include <c4/yml/std/string.hpp>

namespace c4 { namespace yml {
  void read(c4::yml::NodeRef const& n, TH1& h){
//    for(size_t i=0; i<n.num_children(); ++i){
//      T e;
//      n[i]>>e;
//      v->push_back(e);
//    }
  }

    
  void write(c4::yml::NodeRef *bounds, const TAxis& ax){
    if(!ax.IsVariableBinSize()){
      *bounds |= c4::yml::MAP;
      (*bounds)["nbins"] << ax.GetNbins();                
      (*bounds)["min"] << ax.GetXmin();
      (*bounds)["max"] << ax.GetXmax();
    } else {
      *bounds |= c4::yml::SEQ;              
      for(size_t i=1; i<=ax.GetNbins(); ++i){
        bounds->append_child() << ax.GetBinLowEdge(i);      
      }
    }
  }
    
  void write(c4::yml::NodeRef *n, const TH1& h){
    *n |= c4::yml::MAP;
    auto bounds = (*n)["binning"];
    bounds |= c4::yml::MAP;
    auto weights = (*n)["counts"];
    weights |= c4::yml::SEQ;    
    auto errors = (*n)["errors"];    
    errors |= c4::yml::SEQ;    
    if(h.GetDimension()==1){
      write(&bounds,*(h.GetXaxis()));
      for(size_t i=1; i<=h.GetNbinsX(); ++i){
        weights.append_child() << h.GetBinContent(i);
        errors.append_child() << h.GetBinError(i);        
      }    
    } else {
      auto x = bounds["x"];
      write(&x,*(h.GetXaxis()));
      auto y = bounds["y"];
      write(&y,*(h.GetYaxis()));      
      if(h.GetDimension()>2){
        auto z = bounds["z"];
        write(&z,*(h.GetZaxis()));              
      }
      for(size_t i=1; i<=h.GetNbinsX(); ++i){
        auto binx = weights.append_child();
        binx |= c4::yml::SEQ;
        auto binxe = errors.append_child();
        binxe |= c4::yml::SEQ;                                
        for(size_t j=1; j<=h.GetNbinsY(); ++j){
          if(h.GetDimension()>2){      
            auto biny = binx.append_child();
            biny |= c4::yml::SEQ;
            auto binye = binxe.append_child();
            binye |= c4::yml::SEQ;                            
            for(size_t k=1; k<=h.GetNbinsY(); ++k){        
              biny.append_child() << h.GetBinContent(i,j,k);
              binye.append_child() << h.GetBinError(i,j,k);              
            }
          } else {
            binx.append_child() << h.GetBinContent(i,j);
            binxe.append_child() << h.GetBinError(i,j);            
          }
        }
      }          
    }
  }
  }
}
namespace RooStats { namespace HistFactory {
    template<> void RooStats::HistFactory::Sample::Export(c4::yml::NodeRef& n) const {
      auto s = n[c4::to_csubstr(fName)];
      s |= c4::yml::MAP;
      s["type"] << "histogram";
      
      if(fOverallSysList.size() > 0){
        auto overallSys = s["overallSystematics"];
        overallSys |= c4::yml::MAP;
        for(const auto& sys:fOverallSysList){
          auto node = overallSys[c4::to_csubstr(sys.GetName())];
          node |= c4::yml::MAP;        
          node["parameter"] << std::string("alpha_")+sys.GetName();
          node["low"] << sys.GetLow();
          node["high"] << sys.GetHigh();
        }
      }

      if(fNormFactorList.size()>0){
        auto normFactors = s["normFactors"];
        normFactors |= c4::yml::SEQ;
        for(auto& sys:fNormFactorList){
          normFactors.append_child() << sys.GetName();
        }
      }

      if(fHistoSysList.size()>0){
        auto histoSys = s["histogramSystematics"];
        histoSys |= c4::yml::MAP;
        for(size_t i=0; i<fHistoSysList.size(); ++i){
          auto sys = fHistoSysList[i];
          auto node = histoSys[c4::to_csubstr(sys.GetName())];
          node |= c4::yml::MAP;        
          node["parameter"] << std::string("alpha_")+sys.GetName();
          node["dataLow"] << *(sys.GetHistoLow());
          node["dataHigh"] << *(sys.GetHistoHigh());
        }
      }

      // std::vector< RooStats::HistFactory::ShapeSys >    fShapeSysList;
      // std::vector< RooStats::HistFactory::ShapeFactor > fShapeFactorList;
  
      auto tags = s["dict"];
      tags |= c4::yml::MAP;      
      tags["normalizeByTheory"] << fNormalizeByTheory;
      tags["statErrorActivate"] << fStatErrorActivate;

      s["data"] << *(fhNominal.GetObject());
    }
  }
}
#endif


// Some helper functions
// (Not strictly necessary because
//  methods are publicly accessable)


void RooStats::HistFactory::Sample::ActivateStatError() {
  
  fStatError.Activate( true );
  fStatError.SetUseHisto( false );

}


void RooStats::HistFactory::Sample::ActivateStatError( std::string StatHistoName, std::string StatInputFile, std::string StatHistoPath ) {


  fStatError.Activate( true );
  fStatError.SetUseHisto( true );

  fStatError.SetInputFile( StatInputFile );
  fStatError.SetHistoName( StatHistoName );
  fStatError.SetHistoPath( StatHistoPath );

}


void RooStats::HistFactory::Sample::AddOverallSys( std::string SysName, Double_t SysLow, Double_t SysHigh ) {

  RooStats::HistFactory::OverallSys sys;
  sys.SetName( SysName );
  sys.SetLow( SysLow );
  sys.SetHigh( SysHigh );

  fOverallSysList.push_back( sys );

}

void RooStats::HistFactory::Sample::AddOverallSys( const OverallSys& Sys ) {
  fOverallSysList.push_back(Sys);
}

void RooStats::HistFactory::Sample::AddNormFactor( std::string SysName, Double_t SysVal, Double_t SysLow, Double_t SysHigh, bool SysConst ) {

  RooStats::HistFactory::NormFactor norm;

  norm.SetName( SysName );
  norm.SetVal( SysVal );
  norm.SetLow( SysLow );
  norm.SetHigh( SysHigh );
  norm.SetConst( SysConst );

  fNormFactorList.push_back( norm );

}

void RooStats::HistFactory::Sample::AddNormFactor( const NormFactor& Factor ) {
  fNormFactorList.push_back( Factor );
}


void RooStats::HistFactory::Sample::AddHistoSys( std::string SysName, 
std::string SysHistoNameLow,  std::string SysHistoFileLow,  std::string SysHistoPathLow,
						 std::string SysHistoNameHigh, std::string SysHistoFileHigh, std::string SysHistoPathHigh ) {

  RooStats::HistFactory::HistoSys sys;
  sys.SetName( SysName );
  
  sys.SetHistoNameLow( SysHistoNameLow );
  sys.SetHistoPathLow( SysHistoPathLow );
  sys.SetInputFileLow( SysHistoFileLow );

  sys.SetHistoNameHigh( SysHistoNameHigh );
  sys.SetHistoPathHigh( SysHistoPathHigh );
  sys.SetInputFileHigh( SysHistoFileHigh );

  fHistoSysList.push_back( sys );

}

void RooStats::HistFactory::Sample::AddHistoSys( const HistoSys& Sys ) {
  fHistoSysList.push_back( Sys );
}


void RooStats::HistFactory::Sample::AddHistoFactor( std::string SysName, std::string SysHistoNameLow,  std::string SysHistoFileLow,  std::string SysHistoPathLow,  
						    std::string SysHistoNameHigh, std::string SysHistoFileHigh, std::string SysHistoPathHigh ) {

  RooStats::HistFactory::HistoFactor factor;
  factor.SetName( SysName );

  factor.SetHistoNameLow( SysHistoNameLow );
  factor.SetHistoPathLow( SysHistoPathLow );
  factor.SetInputFileLow( SysHistoFileLow );

  factor.SetHistoNameHigh( SysHistoNameHigh );
  factor.SetHistoPathHigh( SysHistoPathHigh );
  factor.SetInputFileHigh( SysHistoFileHigh );

  fHistoFactorList.push_back( factor );

}

void RooStats::HistFactory::Sample::AddHistoFactor( const HistoFactor& Factor ) {
  fHistoFactorList.push_back(Factor);
}


void RooStats::HistFactory::Sample::AddShapeFactor( std::string SysName ) {

  RooStats::HistFactory::ShapeFactor factor;
  factor.SetName( SysName );
  fShapeFactorList.push_back( factor );

}


void RooStats::HistFactory::Sample::AddShapeFactor( const ShapeFactor& Factor ) {
  fShapeFactorList.push_back(Factor);
}


void RooStats::HistFactory::Sample::AddShapeSys( std::string SysName, Constraint::Type SysConstraintType, std::string SysHistoName, std::string SysHistoFile, std::string SysHistoPath ) {

  RooStats::HistFactory::ShapeSys sys;
  sys.SetName( SysName );
  sys.SetConstraintType( SysConstraintType );

  sys.SetHistoName( SysHistoName );
  sys.SetHistoPath( SysHistoPath );
  sys.SetInputFile( SysHistoFile );

  fShapeSysList.push_back( sys );

}

void RooStats::HistFactory::Sample::AddShapeSys( const ShapeSys& Sys ) {
  fShapeSysList.push_back(Sys);
}
