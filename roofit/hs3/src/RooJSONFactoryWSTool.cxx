#include <RooJSONFactoryWSTool.h>

#include <iostream>
#include <fstream>
#include <stdexcept>

#include <RooConstVar.h>
#include <RooRealVar.h>
#include <RooAbsCategory.h>
#include <RooRealProxy.h>
#include <RooStats/ModelConfig.h>

#include "TROOT.h"
#include "RConfigure.h"

#ifdef R__HAS_RYML
#include "RYMLParser.h"
typedef TRYMLTree tree_t;
#else
#include "JSONParser.h"
typedef TJSONTree tree_t;
#endif

ClassImp(RooJSONFactoryWSTool)

namespace {
  // helpers for serializing / deserializing binned datasets
  inline void genIndicesHelper(std::vector<std::vector<int> >& combinations, RooArgList& vars, size_t curridx){
    if(curridx == vars.size()){
      std::vector<int> indices(curridx);
      for(size_t i=0; i<curridx; ++i){
        RooRealVar* v = (RooRealVar*)(vars.at(i));
        indices[i] = v->getBinning().binNumber(v->getVal());
      }
      combinations.push_back(indices);
    } else {
      RooRealVar* v = (RooRealVar*)(vars.at(curridx));
      for(int i=0; i<v->numBins(); ++i){      
        v->setVal(v->getBinning().binCenter(i));
        ::genIndicesHelper(combinations,vars,curridx+1);
      }
    }
  }
}

// maps to hold the importers and exporters for runtime lookup
RooJSONFactoryWSTool::ImportMap RooJSONFactoryWSTool::_importers = RooJSONFactoryWSTool::ImportMap();
RooJSONFactoryWSTool::ExportMap RooJSONFactoryWSTool::_exporters = RooJSONFactoryWSTool::ExportMap();

bool RooJSONFactoryWSTool::registerImporter(const std::string& key, const RooJSONFactoryWSTool::Importer* f){
  if(RooJSONFactoryWSTool::_importers.find(key) != RooJSONFactoryWSTool::_importers.end()) return false;
  RooJSONFactoryWSTool::_importers.insert(std::make_pair(key,f));
  return true;
}
bool RooJSONFactoryWSTool::registerExporter(const TClass* key, const RooJSONFactoryWSTool::Exporter* f){
  if(RooJSONFactoryWSTool::_exporters.find(key) != RooJSONFactoryWSTool::_exporters.end()) return false;
  RooJSONFactoryWSTool::_exporters.insert(std::make_pair(key,f));  
  return true;
}

void RooJSONFactoryWSTool::printImporters(){
  for(const auto& x:RooJSONFactoryWSTool::_importers){
    std::cout << x.first << "\t" << typeid(*x.second).name() << std::endl;
  }
}
void RooJSONFactoryWSTool::printExporters(){
  for(const auto& x:RooJSONFactoryWSTool::_exporters){
    std::cout << x.first << "\t" << typeid(*x.second).name() << std::endl;
  }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// helper functions specific to JSON
///////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
  inline std::string genPrefix(const JSONNode& p,bool trailing_underscore){
    std::string prefix;
    if(!p.is_map()) return prefix;
    if(p.has_child("namespaces")){
      for(const auto& ns:p["namespaces"].children()){
        if(prefix.size() > 0) prefix+="_";
        prefix += ns.val();
      }
    }
    if(trailing_underscore && prefix.size()>0) prefix += "_";
    return prefix;
  }

  inline void importAttributes(RooAbsArg* arg, const JSONNode& n){
    if(!n.is_map()) return;
    if(n.has_child("dict") && n["dict"].is_map()){
      for(const auto& attr:n["dict"].children()){
        arg->setStringAttribute(RooJSONFactoryWSTool::name(attr).c_str(),attr.val().c_str());
      }
    }
    if(n.has_child("tags") && n["tags"].is_seq()){
      for(const auto& attr:n["tags"].children()){
        arg->setAttribute(attr.val().c_str());
      }
    }    
  }
  inline bool checkRegularBins(const TAxis& ax){
    double w = ax.GetXmax() - ax.GetXmin();
    double bw = w / ax.GetNbins();
    for(int i=0; i<=ax.GetNbins(); ++i){
      if( fabs(ax.GetBinUpEdge(i) - (ax.GetXmin()+(bw*i))) > w*1e-6 ) return false;
    }
    return true;
  }
  inline void writeAxis(JSONNode& bounds, const TAxis& ax){
    bool regular = (!ax.IsVariableBinSize()) || checkRegularBins(ax);
    if(regular){
      bounds.set_map();
      bounds["nbins"] << ax.GetNbins();                
      bounds["min"] << ax.GetXmin();
      bounds["max"] << ax.GetXmax();
    } else {
      bounds.set_seq();              
      for(int i=0; i<=ax.GetNbins(); ++i){
        bounds.append_child() << ax.GetBinUpEdge(i);      
      }
    }
  }  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// RooWSFactoryTool expression handling
///////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
  struct JSON_Factory_Expression {
    TClass* tclass;
    std::vector<std::string> arguments;
    std::string generate(const JSONNode& p){
      std::string name(RooJSONFactoryWSTool::name(p));
      std::stringstream expression;
      std::string classname(this->tclass->GetName());
      size_t colon = classname.find_last_of(":");
      if(colon < classname.size()){
        expression << classname.substr(colon+1);
      } else {
        expression << classname;
      }
      expression << "::" << name << "(";
      bool first = true;
      for(auto k:this->arguments){
        if(!first) expression << ",";
        first = false;
        if(k == "true"){
          expression << "1";
          continue;
        } else if(k == "false"){
          expression << "0";
          continue;          
        } else if(!p.has_child(k)){
          std::stringstream err;
          err << "factory expression for class '" << this->tclass->GetName() << "', which expects key '" << k << "' missing from input for object '" << name << "', skipping.";
          RooJSONFactoryWSTool::error(err.str().c_str());
        }
        if(p[k].is_seq()){
          expression << "{";
          bool f = true;
          for(const auto& x:p[k].children()){
            if(!f) expression << ",";
            f=false;
            expression << x.val();
          }
          expression << "}";
        } else {
          expression << p[k].val();
        }          
      }
      expression << ")";
      return expression.str();
    }
  };
  std::map<std::string,JSON_Factory_Expression> _pdfFactoryExpressions;
  std::map<std::string,JSON_Factory_Expression> _funcFactoryExpressions;
}

void RooJSONFactoryWSTool::loadFactoryExpressions(const std::string& fname){
  // load a yml file defining the factory expressions
  std::ifstream infile(fname);
  if(!infile.is_open()){
    std::cerr << "unable to read file '" << fname << "'" << std::endl;
    return;
  }
  try {
    std::cout << "before" << std::endl;
    tree_t p(infile);
    std::cout << "after" << std::endl;    
    const JSONNode& n = p.rootnode();
    for(const auto& cl:n.children()){
      std::string key(RooJSONFactoryWSTool::name(cl));
      if(!cl.has_child("class")){
        std::cerr << "error in file '" << fname << "' for entry '" << key << "': 'class' key is required!" << std::endl;
        continue;
      }
      std::string classname(cl["class"].val());    
      TClass* c = TClass::GetClass(classname.c_str());
      if(!c){
        std::cerr << "unable to find class " << classname << ", skipping." << std::endl;
      } else {
        JSON_Factory_Expression ex;
        ex.tclass = c;
        if(!cl.has_child("arguments")){
          std::cerr << "class " << classname << " seems to have no arguments attached, skipping" << std::endl;
          continue;
        }
        for(const auto& arg:cl["arguments"].children()){
          ex.arguments.push_back(arg.val());
        }
        if(c->InheritsFrom(RooAbsPdf::Class())){
          _pdfFactoryExpressions[key] = ex;
        } else if(c->InheritsFrom(RooAbsReal::Class())){
          _funcFactoryExpressions[key] = ex;        
        } else {
          std::cerr << "class " << classname << " seems to not inherit from any suitable class, skipping" << std::endl;
        }
      }
    }
  } catch (const std::exception& ex){
    std::cout << "caught" << std::endl;        
    std::cerr << "unable to load factory expressions: " << ex.what() << std::endl;
  }  
}
void RooJSONFactoryWSTool::clearFactoryExpressions(){
  // clear all factory expressions
  _pdfFactoryExpressions.clear();
  _funcFactoryExpressions.clear();
}
void RooJSONFactoryWSTool::printFactoryExpressions(){
  // print all factory expressions
  for(auto it:_pdfFactoryExpressions){
    std::cout << it.first;
    std::cout << " " << it.second.tclass->GetName();    
    for(auto v:it.second.arguments){
      std::cout << " " << v;
    }
    std::cout << std::endl;
  }
  for(auto it:_funcFactoryExpressions){
    std::cout << it.first;
    std::cout << " " << it.second.tclass->GetName();    
    for(auto v:it.second.arguments){
      std::cout << " " << v;
    }
    std::cout << std::endl;
  }  
}

std::vector<std::vector<int> > RooJSONFactoryWSTool::generateBinIndices(RooArgList& vars){
  std::vector<std::vector<int> > combinations;
  ::genIndicesHelper(combinations,vars,0);
  return combinations;
}

void  RooJSONFactoryWSTool::exportHistogram(const TH1& h, JSONNode& n, const std::vector<std::string>& varnames){
  n.set_map();
  auto& bounds = n["binning"];
  bounds.set_map();
  auto& weights = n["counts"];
  weights.set_seq();    
  auto& errors = n["errors"];    
  errors.set_seq();
  auto& x = bounds[varnames[0]];
  writeAxis(x,*(h.GetXaxis()));
  if(h.GetDimension()>1){
    auto& y = bounds[varnames[1]];
    writeAxis(y,*(h.GetYaxis()));
    if(h.GetDimension()>2){
      auto& z = bounds[varnames[2]];
      writeAxis(z,*(h.GetZaxis()));              
    }
  }
  for(int i=1; i<=h.GetNbinsX(); ++i){
    if(h.GetDimension()==1){
      weights.append_child() << h.GetBinContent(i);
      errors.append_child() << h.GetBinError(i);
    } else {
      for(int j=1; j<=h.GetNbinsY(); ++j){
        if(h.GetDimension()==2){
          weights.append_child() << h.GetBinContent(i,j);
          errors.append_child() << h.GetBinError(i,j);
        } else {
          for(int k=1; k<=h.GetNbinsZ(); ++k){        
            weights.append_child() << h.GetBinContent(i,j,k);
            errors.append_child() << h.GetBinError(i,j,k);              
          }
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// RooProxy-based export handling
///////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
  struct JSON_Export_Keys {
    std::string type;
    std::map<std::string,std::string> proxies;
  };  
  std::map<TClass*,JSON_Export_Keys> _exportKeys;
}
void RooJSONFactoryWSTool::loadExportKeys(const std::string& fname){
  // load a yml file defining the export keys
  std::ifstream infile(fname);
  if(!infile.is_open()){
    std::cerr << "unable to read file '" << fname << "'" << std::endl;
    return;
  }
  try {
    tree_t p(infile);
    const JSONNode& n = p.rootnode();
    for(const auto& cl:n.children()){
      std::string classname(RooJSONFactoryWSTool::name(cl));
      TClass* c = TClass::GetClass(classname.c_str());
      if(!c){
        std::cerr << "unable to find class " << classname << ", skipping." << std::endl;
      } else {
        JSON_Export_Keys ex;
        if(!cl.has_child("type")){
          std::cerr << "class " << classname << "has not type key set, skipping" << std::endl;
          continue;
        }
        if(!cl.has_child("proxies")){
          std::cerr << "class " << classname << "has no proxies identified, skipping" << std::endl;
          continue;
        }      
        ex.type = cl["type"].val();
        for(const auto& k:cl["proxies"].children()){
          std::string key(RooJSONFactoryWSTool::name(k));
          std::string val(k.val());                
          ex.proxies[key] = val;
        }
        _exportKeys[c] = ex;
      }
    }
  } catch (const std::exception& ex){
    std::cerr << "unable to load export keys: " << ex.what() << std::endl;
  }  
}

void RooJSONFactoryWSTool::clearExportKeys(){
  // clear all export keys
  _exportKeys.clear();
}

void RooJSONFactoryWSTool::printExportKeys(){
  // print all export keys
  for(const auto& it:_exportKeys){
    std::cout << it.first->GetName() << ": " << it.second.type;
    for(const auto& kv:it.second.proxies){
      std::cout << " " << kv.first << "=" << kv.second;
    }
    std::cout << std::endl;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// helper namespace
///////////////////////////////////////////////////////////////////////////////////////////////////////

bool RooJSONFactoryWSTool::find(const JSONNode& n, const std::string& elem){
  // find an attribute
  if(n.is_seq()){
    for(const auto& t:n.children()){
      if(t.val() == elem) return true;
    }
    return false;
  } else if(n.is_map()){
    return n.has_child(elem.c_str());
  }
  return false;
}
  
void RooJSONFactoryWSTool::append(JSONNode& n, const std::string& elem){
    // append an attribute
    n.set_seq();          
    if(!find(n,elem)){
      n.append_child() << elem;
    }
  }
  
void RooJSONFactoryWSTool::exportAttributes(const RooAbsArg* arg, JSONNode& n){
    // export all string attributes of an object
    if(arg->stringAttributes().size() > 0){
      auto& dict = n["dict"];
      dict.set_map();      
      for(const auto& it:arg->stringAttributes()){
        dict[it.first] << it.second;
      }
    }
    if(arg->attributes().size() > 0){
      auto& tags = n["tags"];
      tags.set_seq();      
      for(const auto& it:arg->attributes()){
        RooJSONFactoryWSTool::append(tags,it);
      }
    }
  }

void RooJSONFactoryWSTool::exportVariable(const RooAbsReal*v, JSONNode& n) {
    auto& var = n[v->GetName()];
    const RooConstVar* cv  = dynamic_cast<const RooConstVar*>(v);
    const RooRealVar*  rrv = dynamic_cast<const RooRealVar*>(v);    
    var.set_map();  
    if(cv){
      var["value"] << cv->getVal();
      var["const"] << true;
    } else if(rrv){
      var["value"] << rrv->getVal();
      if(rrv->getMin() > -1e30){
        var["min"] << rrv->getMin();
      }
      if(rrv->getMax() < 1e30){
        var["max"] << rrv->getMax();
      }
      if(rrv->isConstant()){
        var["const"] << rrv->isConstant();
      }
    }
    RooJSONFactoryWSTool::exportAttributes(v,var);
  }
  
void RooJSONFactoryWSTool::exportVariables(const RooArgSet& allElems, JSONNode& n) {
    // export a list of RooRealVar objects
    for(auto* arg:allElems){
      RooRealVar* v = dynamic_cast<RooRealVar*>(arg);
      if(!v) continue;
      exportVariable(v,n);
    }  
  }
  
void RooJSONFactoryWSTool::exportObject(const RooAbsArg* func, JSONNode& n){
    if(func->InheritsFrom(RooConstVar::Class()) && strcmp(func->GetName(),TString::Format("%g",((RooConstVar*)func)->getVal()).Data())==0){
      // for RooConstVar, name and value are the same, so we don't need to do anything
      return;
    } else if(func->InheritsFrom(RooAbsCategory::Class())){
      // categories are created by the respective RooSimultaneous, so we're skipping the export here
      return;
    } else if(func->InheritsFrom(RooRealVar::Class()) || func->InheritsFrom(RooConstVar::Class())){
      auto& vars = n["variables"];
      vars.set_map();
      exportVariable(static_cast<const RooAbsReal*>(func),vars);
      return;
    }
    
    JSONNode& container = func->InheritsFrom(RooAbsPdf::Class()) ? n["pdfs"] : n["functions"] ;
    container.set_map();    
    if(container.has_child(func->GetName())) return;
    
    
    TClass* cl = TClass::GetClass(func->ClassName());
    
    auto it = _exporters.find(cl);
    if(it != _exporters.end()){
      if(it->second->autoExportDependants()) RooJSONFactoryWSTool::exportDependants(func,n);
      auto& elem = container[func->GetName()];
      elem.set_map();
      try {
        if(!it->second->exportObject(this,func,elem)){
          std::cerr << "exporter for type " << cl->GetName() << " does not export objects!" << std::endl;          
        }
        RooJSONFactoryWSTool::exportAttributes(func,elem);        
      } catch (const std::exception& ex){
        std::cerr << ex.what() << ". skipping." << std::endl;
        return;
      }
    } else { // generic import using the factory expressions      
      const auto& dict = _exportKeys.find(cl);
      if(dict == _exportKeys.end()){
        std::cerr << "unable to export class '" << cl->GetName() << "' - no export keys available!" << std::endl;
        std::cerr << "there are several possible reasons for this:" << std::endl;
        std::cerr << " 1. " << cl->GetName() << " is a custom class that you or some package you are using added." << std::endl;
        std::cerr << " 2. " << cl->GetName() << " is a ROOT class that nobody ever bothered to write a serialization definition for." << std::endl;
        std::cerr << " 3. something is wrong with your setup, e.g. you might have called RooJSONFactoryWSTool::clearExportKeys() and/or never successfully read a file defining these keys with RooJSONFactoryWSTool::loadExportKeys(filename)" << std::endl;
        std::cerr << "either way, please make sure that:" << std::endl;
        std::cerr << " 3: you are reading a file with export keys - call RooJSONFactoryWSTool::printExportKeys() to see what is available" << std::endl;
        std::cerr << " 2 & 1: you might need to write a serialization definition yourself. check INSERTLINKHERE to see how to do this!" << std::endl;                              
        return;
      }

      size_t nprox = func->numProxies();
      for(size_t i=0; i<nprox; ++i){
        RooAbsProxy* p = func->getProxy(i);
        
        std::string pname(p->name());
        if(pname[0] == '!') pname.erase(0,1);
        
        auto k = dict->second.proxies.find(pname);
        if(k == dict->second.proxies.end()){
          std::cerr << "failed to find key matching proxy '" << pname << "' for type '" << dict->second.type << "', skipping" << std::endl;
          return;
        }
        
        RooListProxy* l = dynamic_cast<RooListProxy*>(p);
        if(l){
          for(auto e:*l){
            exportObject(e,n);            
          }
          auto& elem = container[func->GetName()];
          elem.set_map();
          elem["type"] << dict->second.type;
          auto& items = elem[k->second];
          items.set_seq();
          for(auto e:*l){          
            items.append_child() << e->GetName();
          }
          RooJSONFactoryWSTool::exportAttributes(func,elem);          
        }
        RooRealProxy* r = dynamic_cast<RooRealProxy*>(p);
        if(r){
          exportObject(&(r->arg()),n);
          auto& elem = container[func->GetName()];
          elem.set_map();
          elem["type"] << dict->second.type;
          elem[k->second] << r->arg().GetName();
          exportAttributes(func,elem);
        }
      }
    }
  }

void RooJSONFactoryWSTool::exportFunctions(const RooArgSet& allElems, JSONNode& n){
    // export a list of functions
    // note: this function assumes that all the dependants of these objects have already been exported
    for(auto* arg:allElems){
      RooAbsReal* func = dynamic_cast<RooAbsReal*>(arg);
      if(!func) continue;
      RooJSONFactoryWSTool::exportObject(func,n);
    }
  }

///////////////////////////////////////////////////////////////////////////////////////////////////////
// importing functions
void RooJSONFactoryWSTool::importFunctions(const JSONNode& n) {
  // import a list of RooAbsReal objects
  if(!n.is_map()) return;
  for(const auto& p:n.children()){
    // some preparations: what type of function are we dealing with here?
    std::string name(RooJSONFactoryWSTool::name(p));
    if(name.size() == 0) continue;
    if(this->_workspace->pdf(name.c_str())) continue;    
    if(!p.is_map()) continue; 
    std::string prefix = ::genPrefix(p,true);
    if(prefix.size() > 0) name = prefix+name;    
    if(!p.has_child("type")){
      coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") no type given for '" << name << "', skipping." << std::endl;
      continue;
    }
    std::string functype(p["type"].val());
    this->importDependants(p);    
    // check for specific implementations
    auto it = _importers.find(functype);
    if(it != _importers.end()){
      try {
        if(!it->second->importFunction(this,p)){
          coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") importer for type " << functype << " does not import functions!" << std::endl;          
        }
      } catch (const std::exception& ex){
        coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") " << ex.what() << ". skipping." << std::endl;
      }
    } else { // generic import using the factory expressions
      auto expr = _funcFactoryExpressions.find(functype);
      if(expr != _funcFactoryExpressions.end()){
        std::string expression = expr->second.generate(p);
        if(!this->_workspace->factory(expression.c_str())){
          coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") failed to create " << expr->second.tclass->GetName() << " '" << name << "', skipping. expression was\n"
                                << expression << std::endl;          
        }
      } else {
        coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") no handling for functype '" << functype << "' implemented, skipping." << "\n"
                              << "there are several possible reasons for this:\n"
                              << " 1. " << functype << " is a custom type that is not available in RooFit.\n"
                              << " 2. " << functype << " is a ROOT class that nobody ever bothered to write a deserialization definition for.\n"
                              << " 3. something is wrong with your setup, e.g. you might have called RooJSONFactoryWSTool::clearFactoryExpressions() and/or never successfully read a file defining these expressions with RooJSONFactoryWSTool::loadFactoryExpressions(filename)\n"
                              << "either way, please make sure that:\n"
                              << " 3: you are reading a file with export keys - call RooJSONFactoryWSTool::printFactoryExpressions() to see what is available\n" 
                              << " 2 & 1: you might need to write a serialization definition yourself. check INSERTLINKHERE to see how to do this!" << std::endl ;
          continue;        
      }
    }
    RooAbsReal* func = this->_workspace->function(name.c_str());
    if(!func){
      coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") something went wrong importing function '" << name << "'." << std::endl;      
    } else {
      ::importAttributes(func,p);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// importing pdfs
void RooJSONFactoryWSTool::importPdfs(const JSONNode& n) {
  // import a list of RooAbsPdf objects
  if(!n.is_map()) return;  
  for(const auto& p:n.children()){
    // general preparations: what type of pdf should we build?
    std::string name(RooJSONFactoryWSTool::name(p));
    if(name.size() == 0) continue;
    if(this->_workspace->pdf(name.c_str())) continue;    
    if(!p.is_map()) continue;
    std::string prefix = ::genPrefix(p,true);
    if(prefix.size() > 0) name = prefix+name;
    if(!p.has_child("type")){
      coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") no type given for '" << name << "', skipping." << std::endl;
      continue;
    }
    bool toplevel = false;
    if(p.has_child("tags")){
      toplevel = RooJSONFactoryWSTool::find(p["tags"],"toplevel");
    }
    std::string pdftype(p["type"].val());
    this->importDependants(p);

    // check for specific implementations
    auto it = _importers.find(pdftype);
    if(it != _importers.end()){
      try {
        if(!it->second->importPdf(this,p)){
          coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") importer for type " << pdftype << " does not import pdfs!" << std::endl;          
        }
      } catch (const std::exception& ex){
        coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") " << ex.what() << ". skipping." << std::endl;
      }
    } else { // default implementation using the factory expressions
      auto expr = _pdfFactoryExpressions.find(pdftype);
      if(expr != _pdfFactoryExpressions.end()){
        std::string expression = expr->second.generate(p);
        if(!this->_workspace->factory(expression.c_str())){
          coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") failed to create " << expr->second.tclass->GetName() << " '" << name << "', skipping. expression was\n"
                                << expression << std::endl;
        }
      } else {
        coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") no handling for pdftype '" << pdftype << "' implemented, skipping." << "\n"
                              << "there are several possible reasons for this:\n"
                              << " 1. " << pdftype << " is a custom type that is not available in RooFit.\n"
                              << " 2. " << pdftype << " is a ROOT class that nobody ever bothered to write a deserialization definition for.\n"
                              << " 3. something is wrong with your setup, e.g. you might have called RooJSONFactoryWSTool::clearFactoryExpressions() and/or never successfully read a file defining these expressions with RooJSONFactoryWSTool::loadFactoryExpressions(filename)\n"
                              << "either way, please make sure that:\n"
                              << " 3: you are reading a file with export keys - call RooJSONFactoryWSTool::printFactoryExpressions() to see what is available\n" 
                              << " 2 & 1: you might need to write a serialization definition yourself. check INSERTLINKHERE to see how to do this!" << std::endl;
        continue;
      }
    }
    // post-processing: make sure that the pdf has been created, and attach needed attributes
    RooAbsPdf* pdf = this->_workspace->pdf(name.c_str());
    if(!pdf){
      coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") something went wrong importing pdf '" << name << "'." << std::endl;      
    } else {
      ::importAttributes(pdf,p);
      if(toplevel){
        // if this is a toplevel pdf, also cereate a modelConfig for it
        std::string mcname = name+"_modelConfig";
        RooStats::ModelConfig* mc = new RooStats::ModelConfig(mcname.c_str(),name.c_str());
        this->_workspace->import(*mc);
        RooStats::ModelConfig* inwsmc = dynamic_cast<RooStats::ModelConfig*>(this->_workspace->obj(mcname.c_str()));
        if(inwsmc){
          inwsmc->SetWS(*(this->_workspace));
          inwsmc->SetPdf(*pdf);
          RooArgSet observables;
          for(auto var:this->_workspace->allVars()){
            if(var->getAttribute("observable")){
              observables.add(*var);
            }
          }
          inwsmc->SetObservables(observables);
        } else {
          coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") object '" << mcname << "' in workspace is not of type RooStats::ModelConfig!" << std::endl;        
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// importing variables
void RooJSONFactoryWSTool::importVariables(const JSONNode& n) {
  // import a list of RooRealVar objects
  if(!n.is_map()) return;  
  for(const auto& p:n.children()){
    std::string name(RooJSONFactoryWSTool::name(p));
    if(this->_workspace->var(name.c_str())) continue;
    if(!p.is_map()){
      coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") node '" << name << "' is not a map, skipping." << std::endl;
      continue;
    }
    double val(p.has_child("value") ? p["value"].val_float() : 1.);
    RooRealVar v(name.c_str(),name.c_str(),val);
    if(p.has_child("min"))    v.setMin     (p["min"].val_float());
    if(p.has_child("max"))    v.setMax     (p["max"].val_float());
    if(p.has_child("nbins"))  v.setBins     (p["nbins"].val_int());    
    if(p.has_child("relErr")) v.setError   (v.getVal()*p["relErr"].val_float());
    if(p.has_child("err"))    v.setError   (           p["err"].val_float());
    if(p.has_child("const"))  v.setConstant(p["const"].val_bool());
    else v.setConstant(false);
    ::importAttributes(&v,p);    
    this->_workspace->import(v);
  }  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// export all dependants (servers) of a RooAbsArg
void RooJSONFactoryWSTool::exportDependants(const RooAbsArg* source, JSONNode& n) {
  // export all the servers of a given RooAbsArg
  auto servers(source->servers());
  for(auto s:servers){
    this->exportObject(s,n);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// import all dependants (servers) of a node
void RooJSONFactoryWSTool::importDependants(const JSONNode& n) {
  // import all the dependants of an object
  if(n.has_child("variables")){
    this->importVariables(n["variables"]);
  }
  if(n.has_child("functions")){
    this->importFunctions(n["functions"]);
  }    
  if(n.has_child("pdfs")){
    this->importPdfs(n["pdfs"]);
  }
}

std::string RooJSONFactoryWSTool::name(const JSONNode& n){
  std::stringstream ss;
  if(n.has_key()){
    ss << n.key();
  } else if(n.is_container()){
    if(n.has_child("name")){
      ss << n["name"].val();
    }
  } else {
    ss << n.val();
  }
  return ss.str();    
}

void RooJSONFactoryWSTool::exportAll( JSONNode& n) {
  // export all ModelConfig objects and attached Pdfs
  RooArgSet main;
  for(auto obj:this->_workspace->allGenericObjects()){
    if(obj->InheritsFrom(RooStats::ModelConfig::Class())){
      RooStats::ModelConfig* mc = static_cast<RooStats::ModelConfig*>(obj);
      RooAbsPdf* pdf = mc->GetPdf();
      this->exportDependants(pdf,n);
      if(mc->GetObservables()){
        auto& vars = n["variables"];
        vars.set_map();                
        for(auto obs:*(mc->GetObservables())){
          RooRealVar* v = dynamic_cast<RooRealVar*>(obs);
          if(v){
            auto& var = vars[obs->GetName()];
            var.set_map();                
            var["nbins"] << v->numBins();
            auto& tags = var["tags"];
            tags.set_seq();
            RooJSONFactoryWSTool::append(tags,"observable");
          }
        }
      }
      main.add(*pdf);     
    }
  }
  for(auto obj:this->_workspace->allPdfs()){
    RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(obj);
    if(!pdf) continue;
    if((pdf->getAttribute("toplevel") || pdf->clients().size()==0) && !main.find(*pdf)){
      this->exportDependants(pdf,n);
      main.add(*pdf);
    }
  }
  if(main.size() > 0){
    auto& pdfs = n["pdfs"];
    pdfs.set_map();  
    RooJSONFactoryWSTool::RooJSONFactoryWSTool::exportFunctions(main,n);
    for(auto& pdf:main){
      auto& node = pdfs[pdf->GetName()];
      node.set_map();
      auto& tags = node["tags"];
      RooJSONFactoryWSTool::append(tags,"toplevel");
    }
  } else {
    std::cerr << "no ModelConfig found in workspace and no pdf identified as toplevel by 'toplevel' attribute or an empty client list. nothing exported!" << std::endl;    
  }
}

Bool_t RooJSONFactoryWSTool::exportJSON( std::ostream& os ) {
  // export the workspace in JSON
  tree_t p;
  JSONNode& n = p.rootnode();    
  n.set_map();
  this->exportAll(n);
  n.writeJSON(os);  
  return true;
}
Bool_t RooJSONFactoryWSTool::exportJSON( const char* filename ) {
  // export the workspace in JSON  
  std::ofstream out(filename);
  if(!out.is_open()){
    coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") invalid output file '" << filename << "'." << std::endl;
    return false;
  }      
  return this->exportJSON(out);
}

Bool_t RooJSONFactoryWSTool::exportYML( std::ostream& os ) {
  // export the workspace in YML
  tree_t p;
  JSONNode& n = p.rootnode();    
  n.set_map();
  this->exportAll(n);
  n.writeYML(os);
  return true;
}
Bool_t RooJSONFactoryWSTool::exportYML( const char* filename ) {
  // export the workspace in YML    
  std::ofstream out(filename);
  if(!out.is_open()){
    coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") invalid output file '" << filename << "'." << std::endl;
    return false;
  }    
  return this->exportYML(out);
}

void RooJSONFactoryWSTool::prepare(){
  gROOT->ProcessLine("using namespace RooStats::HistFactory;");
}

Bool_t RooJSONFactoryWSTool::importJSON( std::istream& is ) {
  // import a JSON file to the workspace
  try {
    tree_t p(is);
    JSONNode& n = p.rootnode();  
    this->prepare();
    this->importDependants(n);
  } catch (const std::exception& ex){
    std::cerr << "unable to import JSON: " << ex.what() << std::endl;
    return false;
  } 
  return true;
}
Bool_t RooJSONFactoryWSTool::importJSON( const char* filename ) {
  // import a JSON file to the workspace  
  std::ifstream infile(filename);
  if(!infile.is_open()){
    coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") invalid input file '" << filename << "'." << std::endl;
    return false;
  }
  return this->importJSON(infile);
}

Bool_t RooJSONFactoryWSTool::importYML( std::istream& is ) {
  // import a YML file to the workspace
  try {
    tree_t p(is);
    JSONNode& n = p.rootnode();
    this->prepare();
    this->importDependants(n);
  } catch (const std::exception& ex){
    std::cerr << "unable to import JSON: " << ex.what() << std::endl;
    return false;
  }    
  return true;
}
Bool_t RooJSONFactoryWSTool::importYML( const char* filename ) {
  // import a YML file to the workspace    
  std::ifstream infile(filename);
  if(!infile.is_open()){
    coutE(InputArguments) << "RooJSONFactoryWSTool(" << GetName() << ") invalid input file '" << filename << "'." << std::endl;
    return false;
  }  
  return this->importYML(infile);
}

RooJSONFactoryWSTool::RooJSONFactoryWSTool(RooWorkspace& ws) : _workspace(&ws){
  // default constructor
}

