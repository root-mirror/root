/**
 *  \ingroup HistFactory 
 */

// A set of utils for navegating HistFactory models
#include <stdexcept>    
#include <typeinfo>

#include "RooStats/HistFactory/ParamHistFunc.h"
#include "TIterator.h"
#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooRealVar.h"
#include "RooProdPdf.h"
#include "TH1.h"

#include "RooStats/HistFactory/HistFactorySimultaneous.h"
#include "RooStats/HistFactory/HistFactoryModelUtils.h"

namespace RooStats{
namespace HistFactory{


  std::string channelNameFromPdf( RooAbsPdf* channelPdf ) {
    std::string channelPdfName = channelPdf->GetName();
    std::string ChannelName = channelPdfName.substr(6, channelPdfName.size() );
    return ChannelName;
  }

  RooAbsPdf* getSumPdfFromChannel( RooAbsPdf* sim_channel ) {

    bool verbose=false;

    if(verbose) std::cout << "Getting the RooRealSumPdf for the channel: " 
			  << sim_channel->GetName() << std::endl;

    std::string channelPdfName = sim_channel->GetName();
    std::string ChannelName = channelPdfName.substr(6, channelPdfName.size() );

    // Now, get the RooRealSumPdf
    // ie the channel WITHOUT constraints
    std::string realSumPdfName = ChannelName + "_model";

    RooAbsPdf* sum_pdf = NULL;        
    TIterator* iter_sum_pdf = sim_channel->getComponents()->createIterator(); //serverIterator();
    bool FoundSumPdf=false;
    RooAbsArg* sum_pdf_arg=NULL;
    while((sum_pdf_arg=(RooAbsArg*)iter_sum_pdf->Next())) {
      std::string NodeClassName = sum_pdf_arg->ClassName();
      if( NodeClassName == std::string("RooRealSumPdf") ) {
	FoundSumPdf=true;
	sum_pdf = (RooAbsPdf*) sum_pdf_arg;
	break;
      }
    }
    if( ! FoundSumPdf ) {
      if(verbose) {
	std::cout << "Failed to find RooRealSumPdf for channel: " << sim_channel->GetName() << std::endl;
	sim_channel->getComponents()->Print("V");
      }
      sum_pdf=NULL;
      //throw std::runtime_error("Failed to find RooRealSumPdf for channel");
    } 
    else {
      if(verbose) std::cout << "Found RooRealSumPdf: " << sum_pdf->GetName() << std::endl;
    }
    delete iter_sum_pdf;
    iter_sum_pdf = NULL;

    return sum_pdf;

  }


 void FactorizeHistFactoryPdf(const RooArgSet &observables, RooAbsPdf &pdf, RooArgList &obsTerms, RooArgList &constraints) {
   // utility function to factorize constraint terms from a pdf 
   // (from G. Petrucciani)
   const std::type_info & id = typeid(pdf);
   if (id == typeid(RooProdPdf)) {
      RooProdPdf *prod = dynamic_cast<RooProdPdf *>(&pdf);
      RooArgList list(prod->pdfList());
      for (int i = 0, n = list.getSize(); i < n; ++i) {
         RooAbsPdf *pdfi = (RooAbsPdf *) list.at(i);
            FactorizeHistFactoryPdf(observables, *pdfi, obsTerms, constraints);
         }
      } else if (id == typeid(RooSimultaneous) || id == typeid(HistFactorySimultaneous) ) {    //|| id == typeid(RooSimultaneousOpt)) {
         RooSimultaneous *sim  = dynamic_cast<RooSimultaneous *>(&pdf);
         RooAbsCategoryLValue *cat = (RooAbsCategoryLValue *) sim->indexCat().Clone();
         for (int ic = 0, nc = cat->numBins((const char *)0); ic < nc; ++ic) {
            cat->setBin(ic);
            FactorizeHistFactoryPdf(observables, *sim->getPdf(cat->getCurrentLabel()), obsTerms, constraints);
         }
         delete cat;
      } else if (pdf.dependsOn(observables)) {
         if (!obsTerms.contains(pdf)) obsTerms.add(pdf);
      } else {
         if (!constraints.contains(pdf)) constraints.add(pdf);
      }
   }

  /*
  void getChannelsFromModel( RooAbsPdf* model, RooArgSet* channels, RooArgSet* channelsWithConstraints ) {

    // Loop through the model
    // Find all channels

    std::string modelClassName = model->ClassName();
    
    if( modelClassName == std::string("RooSimultaneous") || model->InheritsFrom("RooSimultaneous") ) {

      TIterator* simServerItr = model->serverIterator();

      // Loop through the child nodes of the sim pdf
      // and find the channel nodes
      RooAbsArg* sim_channel_arg = NULL;
      while(( sim_channel = (RooAbsArg*) simServerItr->Next() )) {

	RooAbsPdf* sim_channel = (RooAbsPdf*) sim_channel_arg;

	// Ignore the Channel Cat
	std::string channelPdfName = sim_channel->GetName();
	std::string channelClassName = sim_channel->ClassName();
	if( channelClassName == std::string("RooCategory") ) continue;

	// If we got here, we found a channel.
	// Format is model_<ChannelName>

	std::string ChannelName = channelPdfName.substr(6, channelPdfName.size() );

	// Now, get the RooRealSumPdf
	RooAbsPdf* sum_pdf = getSumPdfFromChannel( sim_channel );
	
	
	/ *
	// Now, get the RooRealSumPdf
	// ie the channel WITHOUT constraints

	std::string realSumPdfName = ChannelName + "_model";

	RooAbsPdf* sum_pdf = NULL;        
	TIterator* iter_sum_pdf = sim_channel->getComponents()->createIterator(); //serverIterator();
	bool FoundSumPdf=false;
	RooAbsArg* sum_pdf_arg=NULL;
	while((sum_pdf_arg=(RooAbsArg*)iter_sum_pdf->Next())) {

	  std::string NodeClassName = sum_pdf_arg->ClassName();
	  if( NodeClassName == std::string("RooRealSumPdf") ) {
	    FoundSumPdf=true;
	    sum_pdf = (RooAbsPdf*) sum_pdf_arg;
	    break;
	  }
	}
	if( ! FoundSumPdf ) {
	  std::cout << "Failed to find RooRealSumPdf for channel: " << sim_channel->GetName() << std::endl;
	  sim_channel->getComponents()->Print("V");
	  throw std::runtime_error("Failed to find RooRealSumPdf for channel");
	}                                                  
	delete iter_sum_pdf;
	iter_sum_pdf = NULL;
	* /

	// Okay, now add to the arg sets
	channels->add( *sum_pdf );
	channelsWithConstraints->add( *sim_channel );

      }

      delete simServerItr;

    }
    else {
      std::cout << "Model is not a RooSimultaneous or doesn't derive from one." << std::endl;
      std::cout << "HistFactoryModelUtils isn't yet implemented for these pdf's" << std::endl;
    }

  }
  */

  bool getStatUncertaintyFromChannel( RooAbsPdf* channel, ParamHistFunc*& paramfunc, RooArgList* gammaList ) {

    bool verbose=false;

    // Find the servers of this channel
    //TIterator* iter = channel->serverIterator();
    TIterator* iter = channel->getComponents()->createIterator(); //serverIterator();
    bool FoundParamHistFunc=false;
    RooAbsArg* paramfunc_arg = NULL;        
    while(( paramfunc_arg = (RooAbsArg*) iter->Next() )) {
      std::string NodeName = paramfunc_arg->GetName();
      std::string NodeClassName = paramfunc_arg->ClassName();
      if( NodeClassName != std::string("ParamHistFunc") ) continue;
      if( NodeName.find("mc_stat_") != std::string::npos ) {
	FoundParamHistFunc=true;
	paramfunc = (ParamHistFunc*) paramfunc_arg;
	break;
      }
    }
    if( ! FoundParamHistFunc || !paramfunc ) {
      if(verbose) std::cout << "Failed to find ParamHistFunc for channel: " << channel->GetName() << std::endl;
      return false;
    }                                               
    
    delete iter;
    iter = NULL;

    // Now, get the set of gamma's
    gammaList = (RooArgList*) &( paramfunc->paramList());
    if(verbose) gammaList->Print("V");

    return true;

  }
  

  void getDataValuesForObservables( std::map< std::string, std::vector<double> >& ChannelBinDataMap,
				    RooAbsData* data, RooAbsPdf* pdf ) {

    bool verbose=false;

    //std::map< std::string, std::vector<int>  ChannelBinDataMap;

    RooSimultaneous* simPdf = (RooSimultaneous*) pdf;

    // get category label
    RooArgSet* allobs = (RooArgSet*) data->get();
    TIterator* obsIter = allobs->createIterator();
    RooCategory* cat = NULL;
    RooAbsArg* temp = NULL;
    while( (temp=(RooAbsArg*) obsIter->Next())) {
      // use dynamic cast here instead
      if( strcmp(temp->ClassName(),"RooCategory")==0){
	cat = (RooCategory*) temp;
	break;
      }
    }
    if(verbose) {
      if(!cat) std::cout <<"didn't find category"<< std::endl;
      else std::cout <<"found category"<< std::endl;
    }
    delete obsIter;

    // split dataset
    TList* dataByCategory = data->split(*cat);
    if(verbose) dataByCategory->Print();
    // note :
    // RooAbsData* dataForChan = (RooAbsData*) dataByCategory->FindObject("");

    // loop over channels
    RooCategory* channelCat = (RooCategory*) (&simPdf->indexCat());
    for (const auto& nameIdx : *channelCat) {

      // Get pdf associated with state from simpdf
      RooAbsPdf* pdftmp = simPdf->getPdf(nameIdx.first.c_str());

      std::string ChannelName = pdftmp->GetName(); //tt->GetName();
      if(verbose) std::cout << "Getting data for channel: " << ChannelName << std::endl;
      ChannelBinDataMap[ ChannelName ] = std::vector<double>();

      RooAbsData* dataForChan = (RooAbsData*) dataByCategory->FindObject(nameIdx.first.c_str());
      if(verbose) dataForChan->Print();

      // Generate observables defined by the pdf associated with this state
      RooArgSet* obstmp = pdftmp->getObservables(*dataForChan->get()) ;   
      RooRealVar* obs = ((RooRealVar*)obstmp->first());
      if(verbose) obs->Print();

      //double expected = pdftmp->expectedEvents(*obstmp);

      // set value to desired value (this is just an example)
      // double obsVal = obs->getVal();
      // set obs to desired value of observable
      // obs->setVal( obsVal );
      //double fracAtObsValue = pdftmp->getVal(*obstmp);

      // get num events expected in bin for obsVal
      // double nu = expected * fracAtObsValue;
      
      // an easier way to get n
      TH1* histForN = dataForChan->createHistogram("HhstForN",*obs);
      for(int i=1; i<=histForN->GetNbinsX(); ++i){
	double n = histForN->GetBinContent(i);
	if(verbose) std::cout << "n" <<  i << " = " << n  << std::endl;
	ChannelBinDataMap[ ChannelName ].push_back( n ); 
      }
      delete histForN;
    
    } // End Loop Over Categories
    
    return;

  }


  int getStatUncertaintyConstraintTerm( RooArgList* constraints, RooRealVar* gamma_stat, 
					RooAbsReal*& pois_nom, RooRealVar*& tau ) {
    // Given a set of constraint terms, 
    // find the poisson constraint for the 
    // given gamma and return the mean
    // as well as the 'tau' parameter

    bool verbose=false;

    // To get the constraint term, loop over all constraint terms
    // and look for the gamma_stat name as well as '_constraint'
    // std::string constraintTermName = std::string(gamma_stat->GetName()) + "_constraint";
    TIterator* iter_list = constraints->createIterator();
    RooAbsArg* term_constr=NULL;
    bool FoundConstraintTerm=false;
    RooAbsPdf* constraintTerm=NULL;
    while((term_constr=(RooAbsArg*)iter_list->Next())) {
      std::string TermName = term_constr->GetName();
      // std::cout << "Checking if is a constraint term: " << TermName << std::endl;
      
      //if( TermName.find(gamma_stat->GetName())!=string::npos ) {
      if( term_constr->dependsOn( *gamma_stat) ) {
	if( TermName.find("_constraint")!=std::string::npos ) {
	  FoundConstraintTerm=true;
	  constraintTerm = (RooAbsPdf*) term_constr;
	  break;
	}
      }
    }
    if( FoundConstraintTerm==false ) {
      std::cout << "Error: Couldn't find constraint term for parameter: " << gamma_stat->GetName()
		<< " among constraints: " << constraints->GetName() <<  std::endl;
      constraints->Print("V");
      throw std::runtime_error("Failed to find Gamma ConstraintTerm");
      return -1;
    }
    delete iter_list;

    /*
    RooAbsPdf* constraintTerm = (RooAbsPdf*) constraints->find( constraintTermName.c_str() );
    if( constraintTerm == NULL ) {
      std::cout << "Error: Couldn't find constraint term: " << constraintTermName
		<< " for parameter: " << gamma_stat->GetName()
		<< std::endl;
      throw std::runtime_error("Failed to find Gamma ConstraintTerm");
      return -1;
    }
    */

    // Find the "data" of the poisson term
    // This is the nominal value
    bool FoundNomMean=false;
    TIterator* iter_pois = constraintTerm->serverIterator(); //constraint_args
    RooAbsArg* term_pois ;
    while((term_pois=(RooAbsArg*)iter_pois->Next())) {
      std::string serverName = term_pois->GetName();
      //std::cout << "Checking Server: " << serverName << std::endl;
      if( serverName.find("nom_")!=std::string::npos ) {
	FoundNomMean = true;  
	pois_nom = (RooRealVar*) term_pois;
      }
    }
    if( !FoundNomMean || !pois_nom ) {
      std::cout << "Error: Did not find Nominal Pois Mean parameter in gamma constraint term PoissonMean: "
		<< constraintTerm->GetName() << std::endl;
      throw std::runtime_error("Failed to find Nom Pois Mean");
    }
    else {
      if(verbose) std::cout << "Found Poisson 'data' term: " << pois_nom->GetName() << std::endl;
    }
    delete iter_pois;

    // Taking the constraint term (a Poisson), find
    // the "mean" which is the product: gamma*tau
    // Then, from that mean, find tau
    TIterator* iter_constr = constraintTerm->serverIterator(); //constraint_args
    RooAbsArg* pois_mean_arg=NULL;
    bool FoundPoissonMean = false;
    while(( pois_mean_arg = (RooAbsArg*) iter_constr->Next() )) {
      std::string serverName = pois_mean_arg->GetName();
      if( pois_mean_arg->dependsOn( *gamma_stat ) ) {
	FoundPoissonMean=true;
	// pois_mean = (RooAbsReal*) pois_mean_arg;
	break;
      }
    }
    if( !FoundPoissonMean || !pois_mean_arg ) {
      std::cout << "Error: Did not find PoissonMean parameter in gamma constraint term: "
		<< constraintTerm->GetName() << std::endl;
      throw std::runtime_error("Failed to find PoissonMean");
      return -1;
    }
    else {
      if(verbose) std::cout << "Found Poisson 'mean' term: " << pois_mean_arg->GetName() << std::endl;
    }
    delete iter_constr;


    TIterator* iter_product = pois_mean_arg->serverIterator(); //constraint_args
    RooAbsArg* term_in_product ;
    bool FoundTau=false;
    while((term_in_product=(RooAbsArg*)iter_product->Next())) {
      std::string serverName = term_in_product->GetName();
      //std::cout << "Checking Server: " << serverName << std::endl;
      if( serverName.find("_tau")!=std::string::npos ) {
	FoundTau = true;  
	tau = (RooRealVar*) term_in_product;
      }
    }
    if( !FoundTau || !tau ) {
      std::cout << "Error: Did not find Tau parameter in gamma constraint term PoissonMean: "
		<< pois_mean_arg->GetName() << std::endl;
      throw std::runtime_error("Failed to find Tau");
    }
    else {
      if(verbose) std::cout << "Found Poisson 'tau' term: " << tau->GetName() << std::endl;
    }
    delete iter_product;

    return 0;

  }



} // close RooStats namespace
} // close HistFactory namespace
