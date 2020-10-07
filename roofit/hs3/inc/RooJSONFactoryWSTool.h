#ifndef ROOJSONFACTORYWSTOOL_H
#define ROOJSONFACTORYWSTOOL_H
#include "RooWorkspace.h"
#include "TH1.h"
#include "JSONInterface.h"
#include <string>

class RooJSONFactoryWSTool : public TNamed, public RooPrintable {
 public:
  class Importer {
  public:
    virtual bool importPdf(RooJSONFactoryWSTool*, const JSONNode&) const {
      return false;
    }
    virtual bool importFunction(RooJSONFactoryWSTool*, const JSONNode&) const {
      return false;
    }
    virtual ~Importer(){};
  };
  class Exporter {
  public:
    virtual bool autoExportDependants() const { return true; }
    virtual bool exportObject(RooJSONFactoryWSTool*, const RooAbsArg*, JSONNode&) const {
      return false;
    }
    virtual ~Exporter(){};    
  };   

  typedef std::map<const std::string,const Importer*> ImportMap;
  typedef std::map<const TClass*    ,const Exporter*> ExportMap;  
  
 protected:
  RooWorkspace* _workspace;
  static ImportMap _importers;
  static ExportMap _exporters;    
  void prepare();
 public:
  static std::string name(const JSONNode& n);  
  
  RooJSONFactoryWSTool(RooWorkspace& ws);
  RooWorkspace* workspace() { return this->_workspace; }

  static bool registerImporter(const std::string& key, const RooJSONFactoryWSTool::Importer* f);
  static bool registerExporter(const TClass* key, const RooJSONFactoryWSTool::Exporter* f);
  static void printImporters();
  static void printExporters();

  // error handling helpers
  static void error(const char* s){
    throw std::runtime_error(s);
  }
  static void error(const std::string& s){
    throw std::runtime_error(s);
  }
  template<class T> static std::string concat(const T* items, const std::string& sep=",") {
    // Returns a string being the concatenation of strings in input list <items>
    // (names of objects obtained using GetName()) separated by string <sep>.
    bool first = true;
    std::string text;
    
    // iterate over strings in list
    for(auto it:*items){
      if (!first) {
        // insert separator string
        text += sep;
      } else {
        first = false;
      }
      if(!it) text+="NULL";
      else text+=it->GetName();
    }
    return text;
  }
  template<class T> static std::vector<std::string> names(const T* items) {
    // Returns a string being the concatenation of strings in input list <items>
    // (names of objects obtained using GetName()) separated by string <sep>.
    std::vector<std::string> names;
    // iterate over strings in list
    for(auto it:*items){
      if(!it) names.push_back("NULL");
      else names.push_back(it->GetName());
    }
    return names;
  }
  
  static void exportHistogram(const TH1& h, JSONNode& n, const std::vector<std::string>& obsnames);
  static std::vector<std::vector<int> > generateBinIndices(RooArgList& vars);
  
  Bool_t importJSON(const char* filename);
  Bool_t importYML(const char* filename);
  Bool_t importJSON(std::istream& os);
  Bool_t importYML(std::istream& os);  
  Bool_t exportJSON(const char* fileName) ;
  Bool_t exportYML(const char* fileName) ;
  Bool_t exportJSON(std::ostream& os);
  Bool_t exportYML(std::ostream& os);
  
  static void loadFactoryExpressions(const std::string& fname);
  static void clearFactoryExpressions();
  static void printFactoryExpressions();
  static void loadExportKeys(const std::string& fname);
  static void clearExportKeys();
  static void printExportKeys();  

  void importFunctions(const JSONNode& n);
  void importPdfs(const JSONNode& n);
  void importVariables(const JSONNode& n);
  void importDependants(const JSONNode& n);

  bool find(const JSONNode& n, const std::string& elem);
  void append(JSONNode& n, const std::string& elem);
  
  void exportAttributes(const RooAbsArg* arg, JSONNode& n);
  void exportVariable(const RooAbsReal*v, JSONNode& n);
  void exportVariables(const RooArgSet& allElems, JSONNode& n);
  void exportObject(const RooAbsArg* func, JSONNode& n);
  void exportFunctions(const RooArgSet& allElems, JSONNode& n);  

  void exportAll(JSONNode& n);  
  void exportDependants(const RooAbsArg* source, JSONNode& n);

  ClassDef(RooJSONFactoryWSTool,0)
};
#endif
