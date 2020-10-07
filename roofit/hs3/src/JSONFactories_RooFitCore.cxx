#include <RooJSONFactoryWSTool.h>
#include <JSONInterface.h>

#include <RooDataHist.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooHistFunc.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////
// individually implemented importers
///////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
  class RooProdPdfFactory : public RooJSONFactoryWSTool::Importer {
  public:
    virtual bool importPdf(RooJSONFactoryWSTool* tool, const JSONNode& p) const override {
      std::string name(RooJSONFactoryWSTool::name(p));
      RooArgSet factors;
      if(!p.has_child("factors")){
        RooJSONFactoryWSTool::error("no pdfs of '" + name + "'");
      }
      if(!p["factors"].is_seq()){
        RooJSONFactoryWSTool::error("pdfs '" + name + "' are not a list.");
      }      
      for(const auto& comp:p["factors"].children()){
        std::string pdfname(comp.val());
        RooAbsPdf* pdf = tool->workspace()->pdf(pdfname.c_str());
        if(!pdf){
          RooJSONFactoryWSTool::error("unable to obtain component '" + pdfname + "' of '" + name + "'.");
        }
        factors.add(*pdf);
      }
      RooProdPdf prod(name.c_str(),name.c_str(),factors);
      tool->workspace()->import(prod);
      return true;
    }
  };
  bool _rooprodpdffactory = RooJSONFactoryWSTool::registerImporter("pdfprod",new RooProdPdfFactory());

  class RooAddPdfFactory : public RooJSONFactoryWSTool::Importer {
  public:
    virtual bool importPdf(RooJSONFactoryWSTool* tool, const JSONNode& p) const override {
      std::string name(RooJSONFactoryWSTool::name(p));
      RooArgList pdfs;
      RooArgList coefs;      
      if(!p.has_child("summands")){
        RooJSONFactoryWSTool::error("no summands of '" + name + "'");
      }
      if(!p["summands"].is_seq()){
        RooJSONFactoryWSTool::error("summands '" + name + "' are not a list.");
      }      
      if(!p.has_child("coefficients")){
        RooJSONFactoryWSTool::error("no coefficients of '" + name + "'");
      }
      if(!p["coefficients"].is_seq()){
        RooJSONFactoryWSTool::error("coefficients '" + name + "' are not a list.");
      }      
      for(const auto& comp:p["summands"].children()){
        std::string pdfname(comp.val());
        RooAbsPdf* pdf = tool->workspace()->pdf(pdfname.c_str());
        if(!pdf){
          RooJSONFactoryWSTool::error("unable to obtain component '" + pdfname + "' of '" + name + "'.");
        }
        pdfs.add(*pdf);
      }
      for(const auto& comp:p["coefficients"].children()){
        std::string coefname(comp.val());
        RooAbsArg* coef = tool->workspace()->arg(coefname.c_str());
        if(!coef){
          RooJSONFactoryWSTool::error("unable to obtain component '" + coefname + "' of '" + name + "'.");
        }
        coefs.add(*coef);
      }      
      RooAddPdf add(name.c_str(),name.c_str(),pdfs,coefs);
      tool->workspace()->import(add);
      return true;
    }
  };
  bool _rooaddpdffactory = RooJSONFactoryWSTool::registerImporter("pdfsum",new RooAddPdfFactory());  

  
  class RooSimultaneousFactory : public RooJSONFactoryWSTool::Importer {
  public:
    virtual bool importPdf(RooJSONFactoryWSTool* tool, const JSONNode& p) const override {
      std::string name(RooJSONFactoryWSTool::name(p));
      if(!p.has_child("channels")){
        RooJSONFactoryWSTool::error("no channel components of '" + name + "'");
      }
      tool->importPdfs(p["channels"]);
      std::map< std::string, RooAbsPdf *> components;
      RooCategory cat("channelCat","channelCat");
      for(const auto& comp:p["channels"].children()){
        std::string catname(RooJSONFactoryWSTool::name(comp));
        std::string pdfname(comp.has_val() ? comp.val() : RooJSONFactoryWSTool::name(comp));
        RooAbsPdf* pdf = tool->workspace()->pdf(pdfname.c_str());
        if(!pdf){
          RooJSONFactoryWSTool::error("unable to obtain component '" + pdfname + "' of '" + name + "'");
        }
        components[catname] = pdf;
        cat.defineType(catname.c_str());
      }
      RooSimultaneous simpdf(name.c_str(),name.c_str(),components,cat);
      tool->workspace()->import(simpdf);
      return true;
    }
  };
  bool _roosimultaneousfactory = RooJSONFactoryWSTool::registerImporter("simultaneous",new RooSimultaneousFactory()); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// specialized exporter implementations
///////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
  class RooSimultaneousStreamer : public RooJSONFactoryWSTool::Exporter {
  public:
    bool autoExportDependants() const override { return false; }    
    virtual bool exportObject(RooJSONFactoryWSTool* tool, const RooAbsArg* func, JSONNode& elem) const override {
      const RooSimultaneous* sim = static_cast<const RooSimultaneous*>(func);
      elem["type"] << "simultaneous";
      auto& channels = elem["channels"];
      channels.set_map();
      const auto& indexCat = sim->indexCat();
      for(const auto& cat:indexCat){
        RooAbsPdf* pdf = sim->getPdf(cat->GetName());
        if(!pdf) RooJSONFactoryWSTool::error("no pdf found for category");
        auto& ch = channels[cat->GetName()];
        tool->exportObject(pdf,ch);
      }
      return true;
    }
  };
  bool _roosimultaneousstreamer = RooJSONFactoryWSTool::registerExporter(RooSimultaneous::Class(),new RooSimultaneousStreamer());

  class RooHistFuncStreamer : public RooJSONFactoryWSTool::Exporter {
  public:
    virtual bool exportObject(RooJSONFactoryWSTool*, const RooAbsArg* func, JSONNode& elem) const override {
      const RooHistFunc* hf = static_cast<const RooHistFunc*>(func);
      const RooDataHist& dh = hf->dataHist();
      elem["type"] << "histogram";
      RooArgList vars(*dh.get());
      TH1* hist = hf->createHistogram(RooJSONFactoryWSTool::concat(&vars).c_str());
      auto& data = elem["data"];
      RooJSONFactoryWSTool::exportHistogram(*hist,data,RooJSONFactoryWSTool::names(&vars));
      delete hist;
      return true;
    }
  };
  bool _roohistfuncstreamer = RooJSONFactoryWSTool::registerExporter(RooHistFunc::Class(),new RooHistFuncStreamer());
}


