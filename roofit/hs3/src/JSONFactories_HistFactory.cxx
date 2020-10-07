#include <RooJSONFactoryWSTool.h>

#include <RooStats/HistFactory/ParamHistFunc.h>
#include <RooStats/HistFactory/PiecewiseInterpolation.h>
#include <RooStats/HistFactory/FlexibleInterpVar.h>
#include <RooConstVar.h>
#include <RooCategory.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistFunc.h>
#include <RooRealSumPdf.h>
#include <RooProdPdf.h>
#include <RooSimultaneous.h>
#include <RooPoisson.h>
#include <RooProduct.h>

namespace {
  inline void collectNames(const JSONNode& n,std::vector<std::string>& names){
    for(const auto& c:n.children()){
      names.push_back(RooJSONFactoryWSTool::name(c));
    }
  }

  struct Var {
    int nbins;
    double min;
    double max;
    std::vector<double> bounds;
    
    Var(int n): nbins(n), min(0), max(n) {}
    Var(const JSONNode& val){
      if(val.is_map()){
        if(!val.has_child("nbins")) RooJSONFactoryWSTool::error("no nbins given");
        if(!val.has_child("min"))   RooJSONFactoryWSTool::error("no min given");
        if(!val.has_child("max"))   RooJSONFactoryWSTool::error("no max given");
        this->nbins = val["nbins"].val_int();      
        this->min   = val["min"].val_float();
        this->max   = val["max"].val_float();
      } else if(val.is_seq()){
        for(size_t i=0; i<val.num_children(); ++i){
          this->bounds.push_back(val[i].val_float());
        }
        this->nbins = this->bounds.size();
        this->min = this->bounds[0];
        this->max = this->bounds[this->nbins-1];
      }
    }
  };


    
  
  inline std::map<std::string,Var> readVars(const JSONNode& n,const std::string& obsnamecomp){
    std::map<std::string,Var> vars;
    if(!n.is_map()) return vars;
    if(n.has_child("binning")){
      auto& bounds = n["binning"];
      if(!bounds.is_map()) return vars;    
      if(bounds.has_child("nbins")){
        vars.emplace(std::make_pair("obs_x_"+obsnamecomp,Var(bounds)));
      } else {
        for(const auto& p:bounds.children()){
          vars.emplace(std::make_pair(RooJSONFactoryWSTool::name(p),Var(p)));      
        }
      }
    } else {
      vars.emplace(std::make_pair("obs_x_"+obsnamecomp,Var(n["counts"].num_children())));
    }
    return vars;
  }  
  
  inline void collectObsNames(const JSONNode& n,std::vector<std::string>& obsnames,const std::string& obsnamecomp){
    auto vars = ::readVars(n,obsnamecomp);
    if(obsnames.size() == 0){
      for(auto it:vars){
        obsnames.push_back(it.first);
      }
    }
    if(vars.size() != obsnames.size()){
      throw "inconsistent number of variabels";
    }
  }

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
  
  inline void stackError(const JSONNode& n,std::vector<double>& sumW,std::vector<double>& sumW2){
    if(!n.is_map()) return;    
    if(!n.has_child("counts")) throw "no counts given";
    if(!n.has_child("errors")) throw "no errors given";    
    if(n["counts"].num_children() != n["errors"].num_children()){
      throw "inconsistent bin numbers";
    }
    const size_t nbins = n["counts"].num_children();
    for(size_t ibin=0; ibin<nbins; ++ibin){
      double w = n["counts"][ibin].val_float();
      double e = n["errors"][ibin].val_float();
      if(ibin<sumW.size()) sumW[ibin] += w;
      else sumW.push_back(w);
      if(ibin<sumW2.size()) sumW2[ibin] += e*e;
      else sumW2.push_back(e*e);
    }
  }

  inline RooDataHist* readData(RooWorkspace* ws, const JSONNode& n,const std::string& namecomp,const std::string& obsnamecomp){
    if(!n.is_map()) throw "data is not a map!";
    auto vars = readVars(n,obsnamecomp);
    RooArgList varlist;
    for(auto v:vars){
      std::string name(v.first);
      if(ws->var(name.c_str())){
        varlist.add(*(ws->var(name.c_str())));
      } else {
        auto& var = v.second;
        RooRealVar* rrv = new RooRealVar(name.c_str(),name.c_str(),var.min);
        rrv->setMin(var.min);
        rrv->setMax(var.max);
        rrv->setConstant(true);
        rrv->setBins(var.nbins);
        rrv->setAttribute("observable");
        varlist.addOwned(*rrv);
      }
    }
    RooDataHist* dh = new RooDataHist(("dataHist_"+namecomp).c_str(),namecomp.c_str(),varlist);
    auto bins = RooJSONFactoryWSTool::generateBinIndices(varlist);
    if(!n.has_child("counts")) RooJSONFactoryWSTool::error("no counts given");
    auto& counts = n["counts"];
    if(counts.num_children() != bins.size()) RooJSONFactoryWSTool::error(TString::Format("inconsistent bin numbers: counts=%d, bins=%d",(int)counts.num_children(),(int)(bins.size())));
    for(size_t ibin=0; ibin<bins.size(); ++ibin){
      for(size_t i = 0; i<bins[ibin].size(); ++i){
        RooRealVar* v = (RooRealVar*)(varlist.at(i));
        v->setVal(v->getBinning().binCenter(bins[ibin][i]));
      }
      dh->add(varlist,counts[ibin].val_float());
    }
    return dh;
  }

  RooRealVar* getNP(RooJSONFactoryWSTool* tool, const char* parname) {
    RooRealVar* par = tool->workspace()->var(parname);
    if(!par){
      tool->workspace()->factory(TString::Format("%s[0.,-5,5]",parname).Data());
      par = tool->workspace()->var(parname);
    }
    if(!par) RooJSONFactoryWSTool::error(TString::Format("unable to find nuisance parameter '%s'",parname));
    return par;
  }
  RooAbsPdf* getConstraint(RooJSONFactoryWSTool* tool, const char* sysname) {
    RooAbsPdf* pdf = tool->workspace()->pdf(sysname);
    if(!pdf){
      tool->workspace()->factory(TString::Format("RooGaussian::%s(alpha_%s,0.,1.)",sysname,sysname).Data());
      pdf = tool->workspace()->pdf(sysname);
    }
    if(!pdf) RooJSONFactoryWSTool::error(TString::Format("unable to find constraint term '%s'",sysname));
    return pdf;
  }
  
  class RooHistogramFactory : public RooJSONFactoryWSTool::Importer {
  public:
    virtual bool importFunction(RooJSONFactoryWSTool* tool, const JSONNode& p) const override {
      std::string name(RooJSONFactoryWSTool::name(p));
      std::string prefix = ::genPrefix(p,true);
      if(prefix.size() > 0) name = prefix+name;    
      if(!p.has_child("data")){
        RooJSONFactoryWSTool::error("function '" + name + "' is of histogram type, but does not define a 'data' key");
      }
      try {
        RooArgSet prodElems;
        RooDataHist* dh = ::readData(tool->workspace(),p["data"],name,prefix);       
        RooHistFunc* hf = new RooHistFunc(name.c_str(),RooJSONFactoryWSTool::name(p).c_str(),*(dh->get()),*dh);          
        if(p.has_child("normfactors")){
          for(const auto& nf:p["normfactors"].children()){
            std::string nfname(RooJSONFactoryWSTool::name(nf));
            RooAbsReal* r = tool->workspace()->var(nfname.c_str());
            if(!r){
              RooJSONFactoryWSTool::error("unable to find normalization factor '" + nfname + "'");
            } else {
              prodElems.add(*r);
            }
          }
        }
        if(p.has_child("overallSystematics")){
          RooArgList nps;
          std::vector<double> low;
          std::vector<double> high;
          for(const auto& sys:p["overallSystematics"].children()){
            std::string sysname(RooJSONFactoryWSTool::name(sys));
            std::string parname( sys.has_child("parameter") ? RooJSONFactoryWSTool::name(sys["parameter"]) : "alpha_"+sysname);
            RooAbsReal* par = ::getNP(tool,parname.c_str());
            nps.add(*par);
            low.push_back(sys["low"].val_float());
            high.push_back(sys["high"].val_float());
          }
          RooStats::HistFactory::FlexibleInterpVar* v = new RooStats::HistFactory::FlexibleInterpVar(("overallSys_"+name).c_str(),("overallSys_"+name).c_str(),nps,1.,low,high);
          prodElems.addOwned(*v);
        }
        if(p.has_child("histogramSystematics")){
          RooArgList nps;
          RooArgList low;
          RooArgList high;            
          for(const auto& sys:p["histogramSystematics"].children()){
            std::string sysname(RooJSONFactoryWSTool::name(sys));
            std::string parname( sys.has_child("parameter") ? RooJSONFactoryWSTool::name(sys["parameter"]) : "alpha_"+sysname);            
            RooAbsReal* par = ::getNP(tool,parname.c_str());
            nps.add(*par);
            RooDataHist* dh_low = ::readData(tool->workspace(),p["dataLow"],sysname+"Low_"+name,prefix);
            RooHistFunc hf_low((sysname+"Low_"+name).c_str(),RooJSONFactoryWSTool::name(p).c_str(),*(dh_low->get()),*dh_low);              
            low.add(hf_low);
            RooDataHist* dh_high = ::readData(tool->workspace(),p["dataHigh"],sysname+"High_"+name,prefix);
            RooHistFunc hf_high((sysname+"High_"+name).c_str(),RooJSONFactoryWSTool::name(p).c_str(),*(dh_high->get()),*dh_high);              
            high.add(hf_high);              
          }
          PiecewiseInterpolation* v = new PiecewiseInterpolation(("histoSys_"+name).c_str(),("histoSys_"+name).c_str(),*hf,nps,low,high,true);
          prodElems.addOwned(*v);
        }
        if(prodElems.size() > 0){
          hf->SetName(("hist_"+name).c_str());
          prodElems.addOwned(*hf);          
          RooProduct prod(name.c_str(),name.c_str(),prodElems);
          tool->workspace()->import(prod);
        } else {
          tool->workspace()->import(*hf);
        }
      } catch (const std::runtime_error& e){
        RooJSONFactoryWSTool::error("function '" + name + "' is of histogram type, but 'data' is not a valid definition. " + e.what() + ".");
      }
      return true;
    }
  };
  bool _roohistogramfactory = RooJSONFactoryWSTool::registerImporter("histogram",new RooHistogramFactory());

  class RooRealSumPdfFactory : public RooJSONFactoryWSTool::Importer {
  public:
    virtual bool importPdf(RooJSONFactoryWSTool* tool, const JSONNode& p) const override {
      std::string name(RooJSONFactoryWSTool::name(p));
      RooArgList funcs;
      RooArgList coefs;
      if(!p.has_child("samples")){
        RooJSONFactoryWSTool::error("no samples in '" + name + "', skipping.");
      }
      tool->importFunctions(p["samples"]);      
      RooArgList constraints;
      RooArgList nps;      
      RooConstVar* c = new RooConstVar("1","1",1.);
      std::vector<std::string> usesStatError;
      double statErrorThreshold = 0;
      if(p.has_child("statError")){
        auto& staterr = p["statError"];
        if(staterr.has_child("relThreshold")) statErrorThreshold = staterr["relThreshold"].val_float();
        std::vector<double> relValues;
        if(staterr.has_child("stack")){
          for(const auto& comp:staterr["stack"].children()){
            std::string elem = RooJSONFactoryWSTool::name(comp);
            usesStatError.push_back(elem);
          }
        }
      }
      std::vector<double> sumW;
      std::vector<double> sumW2;
      std::vector<std::string> obsnames;
      std::vector<std::string> sysnames;            
      for(const auto& comp:p["samples"].children()){
        std::string fprefix;
        std::string fname(RooJSONFactoryWSTool::name(comp));
        auto& def = comp.is_container() ? comp : p["functions"][fname.c_str()];
        fprefix = ::genPrefix(def,true);
        if(def["type"].val() == "histogram"){
          try {
            ::collectObsNames(def["data"],obsnames,name);
            if(def.has_child("overallSystematics"))  ::collectNames(def["overallSystematics"],sysnames);
            if(def.has_child("histogramSystematics"))::collectNames(def["histogramSystematics"],sysnames);                                
          } catch (const char* s){
            RooJSONFactoryWSTool::error("function '" + name + "' unable to collect observables from function " + fname + ". " + s );
          }
          if(std::find(usesStatError.begin(),usesStatError.end(),fname) != usesStatError.end()){
            try {              
              ::stackError(def["data"],sumW,sumW2);
            } catch (const char* s){                
              RooJSONFactoryWSTool::error("function '" + name + "' unable to sum statError from function " + fname + ". " + s );
            }                
          }
        }
        RooAbsReal* func = tool->workspace()->function((fprefix+fname).c_str());
        if(!func){
          RooJSONFactoryWSTool::error("unable to obtain component '" + fprefix+fname + "' of '" + name + "'");
        }
        funcs.add(*func);
      }
      RooArgList observables;
      for(auto& obsname:obsnames){
        RooRealVar* obs = tool->workspace()->var(obsname.c_str());
        if(!obs){
          RooJSONFactoryWSTool::error("unable to obtain observable '" + obsname + "' of '" + name + "',");
        }
        observables.add(*obs);
      }
      ParamHistFunc* phf = NULL;
      if(usesStatError.size() > 0){
        RooArgList gammas;
        for(size_t i=0; i<sumW.size(); ++i){
          TString gname = TString::Format("gamma_stat_%s_bin_%d",name.c_str(),(int)i);
          TString tname = TString::Format("tau_stat_%s_bin_%d",name.c_str(),(int)i);
          TString prodname = TString::Format("nExp_stat_%s_bin_%d",name.c_str(),(int)i);
          TString poisname = TString::Format("Constraint_stat_%s_bin_%d",name.c_str(),(int)i);                        
          double tauCV = sumW2[i];
          double err = sqrt(sumW2[i])/sumW[i];
          RooRealVar* g = new RooRealVar(gname.Data(),gname.Data(),1.);
          if(err < statErrorThreshold) g->setConstant(true);
          RooRealVar* tau = new RooRealVar(tname.Data(),tname.Data(),tauCV);
          RooArgSet elems;
          elems.add(*g);
          elems.add(*tau);
          g->setError(err);
          RooProduct* prod = new RooProduct(prodname.Data(),prodname.Data(),elems);
          RooPoisson* pois = new RooPoisson(poisname.Data(),poisname.Data(),*prod,*tau);;        
          gammas.add(*g,true);
          constraints.add(*pois,true);
        }
        nps.add(gammas);
        phf = new ParamHistFunc(TString::Format("%s_mcstat",name.c_str()),"staterror",observables,gammas);
        phf->recursiveRedirectServers(observables);
      }
      for(auto& comp:p["samples"].children()){
        std::string fname(RooJSONFactoryWSTool::name(comp));
        if(std::find(usesStatError.begin(),usesStatError.end(),fname) != usesStatError.end()){
          coefs.add(*phf);
        } else {
          coefs.add(*c);
        }
      }
      for(auto& np:nps){
        for(auto client:np->clients()){
          if(client->InheritsFrom(RooAbsPdf::Class()) && !constraints.find(*client)){
            constraints.add(*client);
          }
        }
      }
      for(auto sysname:sysnames){
        RooAbsPdf* pdf = ::getConstraint(tool,sysname.c_str());
        constraints.add(*pdf);
      }
      if(constraints.getSize() == 0){
        RooRealSumPdf sum(name.c_str(),name.c_str(),funcs,coefs);
        tool->workspace()->import(sum);
      } else {
        RooRealSumPdf sum((name+"_model").c_str(),name.c_str(),funcs,coefs);
        constraints.add(sum);
        RooProdPdf prod(name.c_str(),name.c_str(),constraints);
        tool->workspace()->import(prod);        
      }
      return true;
    }
  };

  bool _roorealsumpdffactory = RooJSONFactoryWSTool::registerImporter("histfactory",new RooRealSumPdfFactory());
}

namespace {
  class FlexibleInterpVarStreamer : public RooJSONFactoryWSTool::Exporter {
  public:
    virtual bool exportObject(RooJSONFactoryWSTool*, const RooAbsArg* func, JSONNode& elem) const override {
      const RooStats::HistFactory::FlexibleInterpVar* fip = static_cast<const RooStats::HistFactory::FlexibleInterpVar*>(func);
      elem["type"] << "interpolation0d";        
      auto& vars = elem["vars"];
      vars.set_seq();
      for(const auto& v:fip->variables()){
        vars.append_child() << v->GetName();
      }
      elem["nom"] << fip->nominal();                
      elem["high"] << fip->high();        
      elem["low"] << fip->low();
      return true;
    }
  };
  bool _flexibleinterpvarstreamer = RooJSONFactoryWSTool::registerExporter(RooStats::HistFactory::FlexibleInterpVar::Class(),new FlexibleInterpVarStreamer());
}


namespace {
  class HistFactoryStreamer : public RooJSONFactoryWSTool::Exporter {
  public:
    bool autoExportDependants() const override { return false; }
    void collectElements(RooArgSet& elems, RooProduct* prod) const {
      for(const auto& e:prod->components()){
        if(e->InheritsFrom(RooProduct::Class())){
          collectElements(elems,(RooProduct*)e);
        } else {
          elems.add(*e);
        }
      }
    }
    bool tryExport(const RooProdPdf* prodpdf, JSONNode& elem) const {
      std::string chname(prodpdf->GetName());
      if(chname.find("model_") == 0){
        chname = chname.substr(6);
      }
      std::vector<RooAbsPdf*> pdfs;
      RooRealSumPdf* sumpdf = NULL;
      for(const auto& v:prodpdf->pdfList()){
        if(v->InheritsFrom(RooRealSumPdf::Class())){
          sumpdf = static_cast<RooRealSumPdf*>(v);
        } else {
          pdfs.push_back(static_cast<RooAbsPdf*>(v));
        }
      }
      if(!sumpdf) return false;
      for(const auto& bw:sumpdf->coefList()){
        if(!bw->InheritsFrom(RooAbsReal::Class())) return false;
      }
      for(const auto& sample:sumpdf->funcList()){
        if(!sample->InheritsFrom(RooProduct::Class())) return false;
      }
      // this seems to be ok
      elem["type"] << "histfactory";
      auto& samples = elem["samples"];
      samples.set_map();
      for(const auto& sample:sumpdf->funcList()){
        std::string samplename = sample->GetName();
        if(samplename.find("L_x_") == 0) samplename = samplename.substr(4);
        auto end = samplename.find("_"+chname);
        if(end < samplename.size()) samplename = samplename.substr(0,end);
        auto& s = samples[samplename];
        s.set_map();
        s["type"] << "histogram";        
        RooProduct* prod = dynamic_cast<RooProduct*>(sample);
        RooArgSet elems;
        collectElements(elems,prod);
        for(const auto& e:elems){
          if(e->InheritsFrom(RooRealVar::Class())){
            auto& norms = s["normFactors"];
            norms.set_seq();
            norms.append_child() << e->GetName();
          } else if(e->InheritsFrom(RooHistFunc::Class())){
            auto& data = s["data"];
            const RooHistFunc* hf = static_cast<const RooHistFunc*>(e);
            const RooDataHist& dh = hf->dataHist();            
            RooArgList vars(*dh.get());
            TH1* hist = hf->createHistogram(RooJSONFactoryWSTool::concat(&vars).c_str());
            RooJSONFactoryWSTool::exportHistogram(*hist,data,RooJSONFactoryWSTool::names(&vars));
            delete hist;
          } else if(e->InheritsFrom(RooStats::HistFactory::FlexibleInterpVar::Class())){
            const RooStats::HistFactory::FlexibleInterpVar* fip = static_cast<const RooStats::HistFactory::FlexibleInterpVar*>(e);
            auto& systs = s["overallSystematics"];
            systs.set_map();
            for(size_t i=0; i<fip->variables().size(); ++i){
              std::string sysname(fip->variables().at(i)->GetName());
              if(sysname.find("alpha_") == 0){
                sysname = sysname.substr(6);
              }
              auto& sys = systs[sysname];
              sys.set_map();                          
              sys["low"] << fip->low()[i];
              sys["high"] << fip->high()[i];
            }
          } else if(e->InheritsFrom(PiecewiseInterpolation::Class())){

          }          
        }
      }
      return true;
    }
    
    virtual bool exportObject(RooJSONFactoryWSTool*, const RooAbsArg* p, JSONNode& elem) const override {
      const RooProdPdf* prodpdf = static_cast<const RooProdPdf*>(p);
      if(tryExport(prodpdf,elem)){
        return true;
      }
      elem["type"] << "pdfprod";        
      auto& factors = elem["factors"];
      factors.set_seq();        
      for(const auto& v:prodpdf->pdfList()){
        factors.append_child() << v->GetName();
      }
      return true;
    }
  };
  bool _histfactorystreamer = RooJSONFactoryWSTool::registerExporter(RooProdPdf::Class(),new HistFactoryStreamer());
}



