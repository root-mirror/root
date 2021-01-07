#ifndef HFJSONTOOL_H
#define HFJSONTOOL_H
#include "RooStats/HistFactory/Measurement.h"
#include "TH1.h"
#include "JSONInterface.h"
#include <string>

namespace RooStats { namespace HistFactory {
class JSONTool {
 protected:
  RooStats::HistFactory::Measurement* _measurement;
  static std::vector<std::string> _strcache;  
  
  void Export(const RooStats::HistFactory::Channel& c, JSONNode& t) const;
  void Export(const RooStats::HistFactory::Sample& s, JSONNode& t) const;  
  
 public:
  JSONTool( RooStats::HistFactory::Measurement* );
  
  void PrintJSON( std::ostream& os = std::cout );    
  void PrintJSON( std::string filename );
  void PrintYAML( std::ostream& os = std::cout );    
  void PrintYAML( std::string filename );  
  void Export(JSONNode& t) const;  

};
  }
}
#endif
