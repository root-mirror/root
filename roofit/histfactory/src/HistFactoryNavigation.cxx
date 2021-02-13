
/** \class RooStats::HistFactory::HistFactoryNavigation
 *  \ingroup HistFactory
 */

#include <iomanip>
#include <sstream>

#include "TFile.h"
#include "TRegexp.h"
#include "TMath.h"

#include "RooRealSumPdf.h"
#include "RooProduct.h"
#include "RooMsgService.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/HistFactory/HistFactoryNavigation.h"
#include "RooStats/HistFactory/HistFactoryException.h"


ClassImp(RooStats::HistFactory::HistFactoryNavigation);


namespace RooStats {
  namespace HistFactory {


    // CONSTRUCTOR
    HistFactoryNavigation::HistFactoryNavigation(ModelConfig* mc) 
      : _minBinToPrint(-1), _maxBinToPrint(-1), 
	_label_print_width(20), _bin_print_width(12) {

      if( !mc ) {
	std::cout << "Error: The supplied ModelConfig is NULL " << std::endl;
	throw hf_exc();
      }

      // Save the model pointer
      RooAbsPdf* pdf_in_mc = mc->GetPdf();
      if( !pdf_in_mc ) {
	std::cout << "Error: The pdf found in the ModelConfig: " << mc->GetName()
		  << " is NULL" << std::endl;
	throw hf_exc();
      }

      // Set the PDF member
      fModel = mc->GetPdf();

      // Get the observables
      RooArgSet* observables_in_mc = const_cast<RooArgSet*>(mc->GetObservables());
      if( !observables_in_mc ) {
	std::cout << "Error: Observable set in the ModelConfig: " << mc->GetName()
		  << " is NULL" << std::endl;
	throw hf_exc();
      }
      if( observables_in_mc->getSize() == 0 ) {
	std::cout << "Error: Observable list: " << observables_in_mc->GetName()
		  << " found in ModelConfig: " << mc->GetName()
		  << " has no entries." << std::endl;
	throw hf_exc();
      }

      // Set the observables member
      fObservables = observables_in_mc;
      
      // Initialize the rest of the members
      _GetNodes(fModel, fObservables);

    }


    // CONSTRUCTOR
    HistFactoryNavigation::HistFactoryNavigation(const std::string& FileName,
						 const std::string& WorkspaceName,
						 const std::string& ModelConfigName) :
      _minBinToPrint(-1), _maxBinToPrint(-1), 
      _label_print_width(20), _bin_print_width(12) {
      
      // Open the File
      TFile* file = new TFile(FileName.c_str());
      if( !file ) {
	std::cout << "Error: Failed to open file: " << FileName << std::endl;
	throw hf_exc();
      }

      // Get the workspace
      RooWorkspace* wspace = (RooWorkspace*) file->Get(WorkspaceName.c_str());
      if( !wspace ) {
	std::cout << "Error: Failed to get workspace: " << WorkspaceName
		  << " from file: " << FileName << std::endl;
	throw hf_exc();
      }

      // Get the ModelConfig
      ModelConfig* mc = (ModelConfig*) wspace->obj(ModelConfigName.c_str());
      if( !mc ) {
	std::cout << "Error: Failed to find ModelConfig: " << ModelConfigName
		  << " from workspace: " << WorkspaceName
		  << " in file: " << FileName << std::endl;
	throw hf_exc();
      }

      // Save the model pointer
      RooAbsPdf* pdf_in_mc = mc->GetPdf();
      if( !pdf_in_mc ) {
	std::cout << "Error: The pdf found in the ModelConfig: " << ModelConfigName
		  << " is NULL" << std::endl;
	throw hf_exc();
      }

      // Set the PDF member
      fModel = pdf_in_mc;

      // Get the observables
      RooArgSet* observables_in_mc = const_cast<RooArgSet*>(mc->GetObservables());
      if( !observables_in_mc ) {
	std::cout << "Error: Observable set in the ModelConfig: " << ModelConfigName
		  << " is NULL" << std::endl;
	throw hf_exc();
      }
      if( observables_in_mc->getSize() == 0 ) {
	std::cout << "Error: Observable list: " << observables_in_mc->GetName()
		  << " found in ModelConfig: " << ModelConfigName
		  << " in file: " << FileName
		  << " has no entries." << std::endl;
	throw hf_exc();
      }

      // Set the observables member
      fObservables = observables_in_mc;
      
      // Initialize the rest of the members
      _GetNodes(fModel, fObservables);

    }


    // CONSTRUCTOR
    HistFactoryNavigation::HistFactoryNavigation(RooAbsPdf* model, RooArgSet* observables) :
      _minBinToPrint(-1), _maxBinToPrint(-1), 
      _label_print_width(20), _bin_print_width(12) {

      // Save the model pointer
      if( !model ) {
	std::cout << "Error: The supplied pdf is NULL" << std::endl;
	throw hf_exc();
      }

      // Set the PDF member
      fModel = model;
      fObservables = observables; 

      // Get the observables
      if( !observables ) {
	std::cout << "Error: Supplied Observable set is NULL" << std::endl;
	throw hf_exc();
      }
      if( observables->getSize() == 0 ) {
	std::cout << "Error: Observable list: " << observables->GetName()
		  << " has no entries." << std::endl;
	throw hf_exc();
      }

      // Initialize the rest of the members
      _GetNodes(fModel, fObservables);

    }


    void HistFactoryNavigation::PrintMultiDimHist(TH1* hist, int bin_print_width) {

      // This is how ROOT makes us loop over histograms :(
      int current_bin = 0;
      int num_bins = hist->GetNbinsX()*hist->GetNbinsY()*hist->GetNbinsZ();
      for(int i = 0; i < num_bins; ++i) {
	// Avoid the overflow/underflow
	current_bin++;
	while( hist->IsBinUnderflow(current_bin) ||
	       hist->IsBinOverflow(current_bin) ) {
	  current_bin++;
	}
	// Check that we should print this bin
	if( _minBinToPrint != -1 && i < _minBinToPrint) continue;
	if( _maxBinToPrint != -1 && i > _maxBinToPrint) break;
	std::cout << std::setw(bin_print_width) << hist->GetBinContent(current_bin);
      }
      std::cout << std::endl;

    }



    RooAbsPdf* HistFactoryNavigation::GetChannelPdf(const std::string& channel) {

      std::map< std::string, RooAbsPdf* >::iterator itr;
      itr = fChannelPdfMap.find(channel);
      
      if( itr == fChannelPdfMap.end() ) {
	std::cout << "Warning: Could not find channel: " << channel
		  << " in pdf: " << fModel->GetName() << std::endl;
	return NULL;
      }
      
      RooAbsPdf* pdf = itr->second;
      if( pdf == NULL ) {
	std::cout << "Warning: Pdf associated with channel: " << channel
		  << " is NULL" << std::endl;
	return NULL;
      }

      return pdf;

    }

    void HistFactoryNavigation::PrintState(const std::string& channel) {

      //int label_print_width = 20;
      //int bin_print_width = 12;
      std::cout << std::endl << channel << ":" << std::endl;

      // Get the map of Samples for this channel
      std::map< std::string, RooAbsReal*> SampleFunctionMap = GetSampleFunctionMap(channel);      

      // Set the size of the print width if necessary
      /*
      for( std::map< std::string, RooAbsReal*>::iterator itr = SampleFunctionMap.begin(); 
	   itr != SampleFunctionMap.end(); ++itr) {
	std::string sample_name = itr->first;
	label_print_width = TMath::Max(label_print_width, (int)sample_name.size()+2);
      }
      */

      // Loop over the SampleFunctionMap and print the individual histograms
      // to get the total histogram for the channel
      int num_bins = 0;
      std::map< std::string, RooAbsReal*>::iterator itr = SampleFunctionMap.begin();
      for( ; itr != SampleFunctionMap.end(); ++itr) {

	std::string sample_name = itr->first;
	std::string tmp_name = sample_name + channel + "_pretty_tmp";
	TH1* sample_hist = GetSampleHist(channel, sample_name, tmp_name);
	num_bins = sample_hist->GetNbinsX()*sample_hist->GetNbinsY()*sample_hist->GetNbinsZ();
	std::cout << std::setw(_label_print_width) << sample_name;

	// Print the content of the histogram
	PrintMultiDimHist(sample_hist, _bin_print_width);
	delete sample_hist;

      }

      // Make the line break as a set of "===============" ...
      std::string line_break;
      int high_bin = _maxBinToPrint==-1 ? num_bins : TMath::Min(_maxBinToPrint, (int)num_bins);
      int low_bin = _minBinToPrint==-1 ? 1 : _minBinToPrint;
      int break_length = (high_bin - low_bin + 1) * _bin_print_width;
      break_length += _label_print_width;
      for(int i = 0; i < break_length; ++i) {
	line_break += "=";
      }
      std::cout << line_break << std::endl;

      std::string tmp_name = channel + "_pretty_tmp";
      TH1* channel_hist = GetChannelHist(channel, tmp_name);
      std::cout << std::setw(_label_print_width) << "TOTAL:";

      // Print the Histogram
      PrintMultiDimHist(channel_hist, _bin_print_width);
      delete channel_hist;

      return;

    }


    void HistFactoryNavigation::PrintState() {
      // Loop over channels and print their states, one after another
      for(unsigned int i = 0; i < fChannelNameVec.size(); ++i) {
	PrintState(fChannelNameVec.at(i));
      }
    }

    
    void HistFactoryNavigation::SetPrintWidths(const std::string& channel) {

      // Get the map of Samples for this channel
      std::map< std::string, RooAbsReal*> SampleFunctionMap = GetSampleFunctionMap(channel);      

      // Get the max of the samples
      for( std::map< std::string, RooAbsReal*>::iterator itr = SampleFunctionMap.begin(); 
	   itr != SampleFunctionMap.end(); ++itr) {
	std::string sample_name = itr->first;
	_label_print_width = TMath::Max(_label_print_width, (int)sample_name.size()+2);
      }

      _label_print_width = TMath::Max( _label_print_width, (int)channel.size() + 7); 
    }


    void HistFactoryNavigation::PrintDataSet(RooDataSet* data, 
					     const std::string& channel_to_print) {
      
      // Print the contents of a 'HistFactory' RooDataset
      // These are stored in a somewhat odd way that makes
      // them difficult to inspect for humans.
      // They have the following layout:
      // =====================================================
      // ChannelA      ChannelB     ChannelCat   Weight
      // -----------------------------------------------------
      // bin_1_center   0           ChannelA     bin_1_height
      // bin_2_center   0           ChannelA     bin_2_height
      //      0        bin_1_center ChannelB     bin_1_height
      //      0        bin_2_center ChannelB     bin_2_height
      //                        ...etc...
      // =====================================================

      // int label_print_width = 20;
      // int bin_print_width = 12;

      // Get the Data Histogram for this channel
      for( unsigned int i_chan=0; i_chan < fChannelNameVec.size(); ++i_chan) {

	std::string channel_name = fChannelNameVec.at(i_chan);

	// If we pass a channel string, we only print that one channel
	if( channel_to_print != "" && channel_name != channel_to_print) continue;

	TH1* data_hist = GetDataHist(data, channel_name, channel_name+"_tmp");
	std::cout << std::setw(_label_print_width) << channel_name + " (data)";

	// Print the Histogram
	PrintMultiDimHist(data_hist, _bin_print_width);
	delete data_hist;
      }
    }


    void HistFactoryNavigation::PrintModelAndData(RooDataSet* data) {
      // Loop over all channels and print model
      // (including all samples) and compare
      // it to the supplied dataset

      for( unsigned int i = 0; i < fChannelNameVec.size(); ++i) {
	std::string channel = fChannelNameVec.at(i);
	SetPrintWidths(channel);
	PrintState(channel);
	PrintDataSet(data, channel);
      }
      
      std::cout << std::endl;

    }


    void HistFactoryNavigation::PrintParameters(bool IncludeConstantParams) {

      // Get the list of parameters
      RooArgSet* params = fModel->getParameters(*fObservables);
      
      std::cout << std::endl;

      // Create the title row
      std::cout << std::setw(30) << "Parameter";
      std::cout << std::setw(15) << "Value"
		<< std::setw(15) << "Error Low" 
		<< std::setw(15) << "Error High"
		<< std::endl;
      
      // Loop over the parameters and print their values, etc
      TIterator* paramItr = params->createIterator();
      RooRealVar* param = NULL;
      while( (param=(RooRealVar*)paramItr->Next()) ) {

	if( !IncludeConstantParams && param->isConstant() ) continue;

	std::cout << std::setw(30) << param->GetName();
	std::cout << std::setw(15) << param->getVal();
	if( !param->isConstant() ) {
	  std::cout << std::setw(15) << param->getErrorLo() << std::setw(15) << param->getErrorHi();
	}
	std::cout<< std::endl;
      }
      
      std::cout << std::endl;

      return;
    }

    void HistFactoryNavigation::PrintChannelParameters(const std::string& channel,
						       bool IncludeConstantParams) {

      // Get the list of parameters
      RooArgSet* params = fModel->getParameters(*fObservables);

      // Get the pdf for this channel
      RooAbsPdf* channel_pdf = GetChannelPdf(channel);

      std::cout << std::endl;

      // Create the title row
      std::cout << std::setw(30) << "Parameter";
      std::cout << std::setw(15) << "Value"
		<< std::setw(15) << "Error Low" 
		<< std::setw(15) << "Error High"
		<< std::endl;
      
      // Loop over the parameters and print their values, etc
      TIterator* paramItr = params->createIterator();
      RooRealVar* param = NULL;
      while( (param=(RooRealVar*)paramItr->Next()) ) {

	if( !IncludeConstantParams && param->isConstant() ) continue;

	if( findChild(param->GetName(), channel_pdf)==NULL ) continue;

	std::cout << std::setw(30) << param->GetName();
	std::cout << std::setw(15) << param->getVal();
	if( !param->isConstant() ) {
	  std::cout << std::setw(15) << param->getErrorLo() << std::setw(15) << param->getErrorHi();
	}
	std::cout<< std::endl;
      }
      
      std::cout << std::endl;

      return;
    }


    void HistFactoryNavigation::PrintSampleParameters(const std::string& channel,
						      const std::string& sample,
						      bool IncludeConstantParams) {

      // Get the list of parameters
      RooArgSet* params = fModel->getParameters(*fObservables);

      // Get the pdf for this channel
      RooAbsReal* sample_func = SampleFunction(channel, sample);
      
      std::cout << std::endl;

      // Create the title row
      std::cout << std::setw(30) << "Parameter";
      std::cout << std::setw(15) << "Value"
		<< std::setw(15) << "Error Low" 
		<< std::setw(15) << "Error High"
		<< std::endl;
      
      // Loop over the parameters and print their values, etc
      TIterator* paramItr = params->createIterator();
      RooRealVar* param = NULL;
      while( (param=(RooRealVar*)paramItr->Next()) ) {

	if( !IncludeConstantParams && param->isConstant() ) continue;

	if( findChild(param->GetName(), sample_func)==NULL ) continue;

	std::cout << std::setw(30) << param->GetName();
	std::cout << std::setw(15) << param->getVal();
	if( !param->isConstant() ) {
	  std::cout << std::setw(15) << param->getErrorLo() << std::setw(15) << param->getErrorHi();
	}
	std::cout<< std::endl;
      }
      
      std::cout << std::endl;

      return;
    }



    double HistFactoryNavigation::GetBinValue(int bin, const std::string& channel) {
      // Get the total bin height for the ith bin (ROOT indexing convention)
      // in channel 'channel'
      // (Could be optimized, it uses an intermediate histogram for now...)
      
      // Get the histogram, fetch the bin content, and return
      TH1* channel_hist_tmp = GetChannelHist(channel, (channel+"_tmp").c_str());
      double val = channel_hist_tmp->GetBinContent(bin);
      delete channel_hist_tmp;
      return val;
    }


    double HistFactoryNavigation::GetBinValue(int bin, const std::string& channel, const std::string& sample){  
      // Get the total bin height for the ith bin (ROOT indexing convention)
      // in channel 'channel'
      // (This will be slow if you plan on looping over it.
      //  Could be optimized, it uses an intermediate histogram for now...)

      // Get the histogram, fetch the bin content, and return
      TH1* sample_hist_tmp = GetSampleHist(channel, sample,  (channel+"_tmp").c_str());
      double val = sample_hist_tmp->GetBinContent(bin);
      delete sample_hist_tmp;
      return val;
    }


    std::map< std::string, RooAbsReal*> HistFactoryNavigation::GetSampleFunctionMap(const std::string& channel) {
      // Get a map of strings to function pointers, 
      // which each function cooresponds to a sample

      std::map< std::string, std::map< std::string, RooAbsReal*> >::iterator channel_itr;
      channel_itr = fChannelSampleFunctionMap.find(channel);
      if( channel_itr==fChannelSampleFunctionMap.end() ){
	std::cout << "Error: Channel: " << channel << " not found in Navigation" << std::endl;
	throw hf_exc();
      }

      return channel_itr->second;
    }


    RooAbsReal* HistFactoryNavigation::SampleFunction(const std::string& channel, const std::string& sample){
      // Return the function object pointer cooresponding
      // to a particular sample in a particular channel

      std::map< std::string, std::map< std::string, RooAbsReal*> >::iterator channel_itr;
      channel_itr = fChannelSampleFunctionMap.find(channel);
      if( channel_itr==fChannelSampleFunctionMap.end() ){
	std::cout << "Error: Channel: " << channel << " not found in Navigation" << std::endl;
	throw hf_exc();
      }

      std::map< std::string, RooAbsReal*>& SampleMap = channel_itr->second;
      std::map< std::string, RooAbsReal*>::iterator sample_itr;
      sample_itr = SampleMap.find(sample);
      if( sample_itr==SampleMap.end() ){
	std::cout << "Error: Sample: " << sample << " not found in Navigation" << std::endl;
	throw hf_exc();
      }
      
      return sample_itr->second;

    }


    RooArgSet* HistFactoryNavigation::GetObservableSet(const std::string& channel) {
      // Get the observables for a particular channel

      std::map< std::string, RooArgSet*>::iterator channel_itr;
      channel_itr = fChannelObservMap.find(channel);
      if( channel_itr==fChannelObservMap.end() ){
	std::cout << "Error: Channel: " << channel << " not found in Navigation" << std::endl;
	throw hf_exc();
      }
      
      return channel_itr->second;

    }


    TH1* HistFactoryNavigation::GetSampleHist(const std::string& channel, const std::string& sample,
					      const std::string& hist_name) {
      // Get a histogram of the expected values for
      // a particular sample in a particular channel
      // Give a name, or a default one will be used

      RooArgList observable_list( *GetObservableSet(channel) );
      
      std::string name = hist_name;
      if(hist_name=="") name = channel + "_" + sample + "_hist";

      RooAbsReal* sample_function = SampleFunction(channel, sample);

      return MakeHistFromRooFunction( sample_function, observable_list, name );
				     
    }


    TH1* HistFactoryNavigation::GetChannelHist(const std::string& channel, const std::string& hist_name) {
      // Get a histogram of the total expected value
      // per bin for this channel
      // Give a name, or a default one will be used

      RooArgList observable_list( *GetObservableSet(channel) );
      
      std::map< std::string, RooAbsReal*> SampleFunctionMap = GetSampleFunctionMap(channel);

      // Okay, 'loop' once 
      TH1* total_hist=NULL;
      std::map< std::string, RooAbsReal*>::iterator itr = SampleFunctionMap.begin();
      for( ; itr != SampleFunctionMap.end(); ++itr) {
	std::string sample_name = itr->first;
	std::string tmp_hist_name = sample_name + "_hist_tmp";
	RooAbsReal* sample_function = itr->second;
	TH1* sample_hist = MakeHistFromRooFunction(sample_function, observable_list, 
						   tmp_hist_name);
	total_hist = (TH1*) sample_hist->Clone("TotalHist");
	delete sample_hist;
	break;
      }
      total_hist->Reset();

      // Loop over the SampleFunctionMap and add up all the histograms
      // to get the total histogram for the channel
      itr = SampleFunctionMap.begin();
      for( ; itr != SampleFunctionMap.end(); ++itr) {
	std::string sample_name = itr->first;
	std::string tmp_hist_name = sample_name + "_hist_tmp";
	RooAbsReal* sample_function = itr->second;
	TH1* sample_hist = MakeHistFromRooFunction(sample_function, observable_list, 
						   tmp_hist_name);
	total_hist->Add(sample_hist);
	delete sample_hist;
      }

      if(hist_name=="") total_hist->SetName(hist_name.c_str());
      else total_hist->SetName( (channel + "_hist").c_str() ); 

      return total_hist;

    }


    std::vector< std::string > HistFactoryNavigation::GetChannelSampleList(const std::string& channel) {
      
      std::vector<std::string> sample_list;

      std::map< std::string, RooAbsReal*> sample_map = fChannelSampleFunctionMap[channel];
      std::map< std::string, RooAbsReal*>::iterator itr = sample_map.begin();;
      for( ; itr != sample_map.end(); ++itr) {
	sample_list.push_back( itr->first );
      }

      return sample_list;

    }

    
    THStack* HistFactoryNavigation::GetChannelStack(const std::string& channel, 
						    const std::string& name) {
    
      THStack* stack = new THStack(name.c_str(), "");

      std::vector< std::string > samples = GetChannelSampleList(channel);

      // Add the histograms
      for( unsigned int i=0; i < samples.size(); ++i) {
	std::string sample_name = samples.at(i);
	TH1* hist = GetSampleHist(channel, sample_name, sample_name+"_tmp");
	hist->SetLineColor(2+i);
	hist->SetFillColor(2+i);
	stack->Add(hist);
      }

      return stack;

    }

  
    TH1* HistFactoryNavigation::GetDataHist(RooDataSet* data, const std::string& channel, 
					    const std::string& name) {
					    
      // TO DO:
      // MAINTAIN THE ACTUAL RANGE, USING THE OBSERVABLES
      // MAKE IT WORK FOR MULTI-DIMENSIONAL
      // 

      // If the dataset covers multiple categories,
      // Split the dataset based on the categories
      if(strcmp(fModel->ClassName(),"RooSimultaneous")==0){

	// If so, get a list of the component pdf's:
	RooSimultaneous* simPdf = (RooSimultaneous*) fModel;
	RooCategory* channelCat = (RooCategory*) (&simPdf->indexCat());

	TList* dataset_list = data->split(*channelCat);

	data = dynamic_cast<RooDataSet*>( dataset_list->FindObject(channel.c_str()) );
	
      }

      RooArgList vars( *GetObservableSet(channel) );

      int dim = vars.getSize();

      TH1* hist = NULL;

      if( dim==1 ) {
	RooRealVar* varX = (RooRealVar*) vars.at(0);
	hist = data->createHistogram( name.c_str(),*varX, RooFit::Binning(varX->getBinning()) );
      }
      else if( dim==2 ) {
	RooRealVar* varX = (RooRealVar*) vars.at(0);
	RooRealVar* varY = (RooRealVar*) vars.at(1);
	hist = data->createHistogram( name.c_str(),*varX, RooFit::Binning(varX->getBinning()),
				      RooFit::YVar(*varY, RooFit::Binning(varY->getBinning())) );
      }
      else if( dim==3 ) {
	RooRealVar* varX = (RooRealVar*) vars.at(0);
	RooRealVar* varY = (RooRealVar*) vars.at(1);
	RooRealVar* varZ = (RooRealVar*) vars.at(2);
	hist = data->createHistogram( name.c_str(),*varX, RooFit::Binning(varX->getBinning()), 
				      RooFit::YVar(*varY, RooFit::Binning(varY->getBinning())),
				      RooFit::YVar(*varZ, RooFit::Binning(varZ->getBinning())) );
      }
      else {
	std::cout << "Error: To Create Histogram from RooDataSet, Dimension must be 1, 2, or 3" << std::endl;
	std::cout << "Observables: " << std::endl;
	vars.Print("V");
	throw hf_exc();
      }

      return hist;

    }


    void HistFactoryNavigation::DrawChannel(const std::string& channel, RooDataSet* data) {
    
      // Get the stack
      THStack* stack = GetChannelStack(channel, channel+"_stack_tmp");

      stack->Draw();
      
      if( data!=NULL ) {
	TH1* data_hist = GetDataHist(data, channel, channel+"_data_tmp");
	data_hist->Draw("SAME");
      }

    }
  


    RooArgSet HistFactoryNavigation::_GetAllProducts(RooProduct* node) {

      // An internal method to recursively get all products,
      // including if a RooProduct is a Product of RooProducts
      // etc

      RooArgSet allTerms;

      // Get All Subnodes of this product
      RooArgSet productComponents = node->components();
      
      // Loop over the subnodes and add
      TIterator* itr = productComponents.createIterator();
      RooAbsArg* arg = NULL;
      while( (arg=(RooAbsArg*)itr->Next()) ) {
	std::string ClassName = arg->ClassName();
	if( ClassName == "RooProduct" ) {
	  RooProduct* prod = dynamic_cast<RooProduct*>(arg);
	  allTerms.add( _GetAllProducts(prod) );
	}
	else {
	  allTerms.add(*arg);
	}
      }
      delete itr;

      return allTerms;

    }




    void HistFactoryNavigation::_GetNodes(RooAbsPdf* modelPdf, const RooArgSet* observables) {

      // Get the pdf from the ModelConfig
      //RooAbsPdf* modelPdf = mc->GetPdf();
      //RooArgSet* observables = mc->GetObservables();

      // Create vectors to hold the channel pdf's
      // as well as the set of observables for each channel
      //std::map< std::string, RooAbsPdf* >  channelPdfMap;
      //std::map< std::string, RooArgSet* >  channelObservMap;
      
      // Check if it is a simultaneous pdf or not
      // (if it's an individual channel, it won't be, if it's
      // combined, it's simultaneous)
      // Fill the channel vectors based on the structure
      // (Obviously, if it's not simultaneous, there will be
      // only one entry in the vector for the single channel)
      if(strcmp(modelPdf->ClassName(),"RooSimultaneous")==0){

	// If so, get a list of the component pdf's:
	RooSimultaneous* simPdf = (RooSimultaneous*) modelPdf;
	RooCategory* channelCat = (RooCategory*) (&simPdf->indexCat());

	// Iterate over the categories and get the
	// pdf and observables for each category
	for (const auto& nameIdx : *channelCat) {
	  const std::string& ChannelName = nameIdx.first;
	  fChannelNameVec.push_back( ChannelName );
	  RooAbsPdf* pdftmp = simPdf->getPdf(ChannelName.c_str()) ;
	  RooArgSet* obstmp = pdftmp->getObservables(*observables) ;
	  fChannelPdfMap[ChannelName] = pdftmp;
	  fChannelObservMap[ChannelName] =  obstmp;
	}

      } else { 
	RooArgSet* obstmp = modelPdf->getObservables(*observables) ;	
	// The channel name is model_CHANNEL
	std::string ChannelName = modelPdf->GetName();
	ChannelName = ChannelName.replace(0, 6, "");
	fChannelNameVec.push_back(ChannelName);
	fChannelPdfMap[ChannelName] = modelPdf;
	fChannelObservMap[ChannelName] = obstmp;

      }

      // Okay, now we have maps of the pdfs 
      // and the observable list per channel
      // We then loop over the channel pdfs:
      // and find their RooRealSumPdfs
      // std::map< std::string, RooRealSumPdf* > channelSumNodeMap;

      for( unsigned int i = 0; i < fChannelNameVec.size(); ++i ) {

	std::string ChannelName = fChannelNameVec.at(i);
	RooAbsPdf* pdf = fChannelPdfMap[ChannelName];
	//std::string Name = fChannelNameMap[ChannelName];

	// Loop over the pdf's components and find
	// the (one) that is a RooRealSumPdf
	// Based on the mode, we assume that node is 
	// the "unconstrained" pdf node for that channel
	RooArgSet* components = pdf->getComponents();
	TIterator* argItr = components->createIterator();
	RooAbsArg* arg = NULL;
	while( (arg=(RooAbsArg*)argItr->Next()) ) {
	  std::string ClassName = arg->ClassName();
	  if( ClassName == "RooRealSumPdf" ) {
	    fChannelSumNodeMap[ChannelName] = (RooRealSumPdf*) arg;
	    break;
	  }
	}
      }
      
      // Okay, now we have all necessary
      // nodes filled for each channel.
      for( unsigned int i = 0; i < fChannelNameVec.size(); ++i ) {

	std::string ChannelName = fChannelNameVec.at(i);
	RooRealSumPdf* sumPdf = dynamic_cast<RooRealSumPdf*>(fChannelSumNodeMap[ChannelName]);
	
	// We now take the RooRealSumPdf and loop over
	// its component functions.  The RooRealSumPdf turns
	// a list of functions (expected events or bin heights
	// per sample) and turns it into a pdf.
	// Therefore, we loop over it to find the expected
	// height for the various samples

	// First, create a map to store the function nodes
	// for each sample in this channel
	std::map< std::string, RooAbsReal*> sampleFunctionMap;

	// Loop over the sample nodes in this
	// channel's RooRealSumPdf
	RooArgList nodes = sumPdf->funcList();
	TIterator* sampleItr = nodes.createIterator();
	RooAbsArg* sample;
	while( (sample=(RooAbsArg*)sampleItr->Next()) ) {

	  // Cast this node as a function
	  RooAbsReal* func = (RooAbsReal*) sample;

	  // Do a bit of work to get the name of each sample
	  std::string SampleName = sample->GetName();
	  if( SampleName.find("L_x_") != std::string::npos ) {
	    size_t index = SampleName.find("L_x_");
	    SampleName.replace( index, 4, "" );
	  }
	  if( SampleName.find(ChannelName.c_str()) != std::string::npos ) {
	    size_t index = SampleName.find(ChannelName.c_str());
	    SampleName = SampleName.substr(0, index-1);
	  }

	  // And simply save this node into our map
	  sampleFunctionMap[SampleName] = func;

	}

	fChannelSampleFunctionMap[ChannelName] = sampleFunctionMap;

	// Okay, now we have a list of histograms
	// representing the samples for this channel.

      }

    }


    RooAbsArg* HistFactoryNavigation::findChild(const std::string& name, RooAbsReal* parent) const {
      
      RooAbsArg* term=NULL;

      // Check if it is a "component",
      // ie a sub node:
      RooArgSet* components = parent->getComponents();
      TIterator* argItr = components->createIterator();
      RooAbsArg* arg = NULL;
      while( (arg=(RooAbsArg*)argItr->Next()) ) {
	std::string ArgName = arg->GetName();
	if( ArgName == name ) {
	  term = arg; //dynamic_cast<RooAbsReal*>(arg);
	  break;
	}
      }
      delete components;
      delete argItr;

      if( term != NULL ) return term;

      // If that failed, 
      // Check if it's a Parameter
      // (ie a RooRealVar)
      RooArgSet* args = new RooArgSet();
      RooArgSet* paramSet = parent->getParameters(args);
      TIterator* paramItr = paramSet->createIterator();
      RooAbsArg* param = NULL;
      while( (param=(RooAbsArg*)paramItr->Next()) ) {
	std::string ParamName = param->GetName();
	if( ParamName == name ) {
	  term = param; //dynamic_cast<RooAbsReal*>(arg);
	  break;
	}
      }
      delete args;
      delete paramSet;
      delete paramItr;

      /* Not sure if we want to be silent
	 But since we're returning a pointer which can be NULL,
	 I think it's the user's job to do checks on it.
	 A dereference will always cause a crash, so it won't
	 be silent for long...
	 if( term==NULL ) {
	 std::cout << "Error: Failed to find node: " << name
	 << " as a child of: " << parent->GetName()
	 << std::endl;
	 }
      */

      return term;

    }


    RooAbsReal* HistFactoryNavigation::GetConstraintTerm(const std::string& parameter) {
     
      std::string ConstraintTermName = parameter + "Constraint";

      // First, as a sanity check, let's see if the parameter
      // itself actually exists and if the model depends on it:
      RooRealVar* param = dynamic_cast<RooRealVar*>(findChild(parameter, fModel));
      if( param==NULL ) {
	std::cout << "Error: Couldn't Find parameter: " << parameter << " in model."
		  << std::endl;
	return NULL;
      }

      // The "gamma"'s use a different constraint term name
      if( parameter.find("gamma_stat_") != std::string::npos ) { 
	ConstraintTermName = parameter + "_constraint";
      }

      // Now, get the constraint itself
      RooAbsReal* term = dynamic_cast<RooAbsReal*>(findChild(ConstraintTermName, fModel));

      if( term==NULL ) {
	std::cout << "Error: Couldn't Find constraint term for parameter: " << parameter
		  << " (Looked for '" << ConstraintTermName << "')" << std::endl;
	return NULL;
      }

      return term;

    }


    double HistFactoryNavigation::GetConstraintUncertainty(const std::string& parameter) {
      
      RooAbsReal* constraintTerm = GetConstraintTerm(parameter);
      if( constraintTerm==NULL ) {
	std::cout << "Error: Cannot get uncertainty because parameter: " << parameter
		  << " has no constraint term" << std::endl;
	throw hf_exc();
      }

      // Get the type of constraint
      std::string ConstraintType = constraintTerm->IsA()->GetName();

      // Find its value
      double sigma = 0.0;

      if( ConstraintType == "" ) {
	std::cout << "Error: Constraint type is an empty string."  
		  << " This simply should not be." << std::endl;
	throw hf_exc();
      }
      else if( ConstraintType == "RooGaussian" ){
	
	// Gaussian errors are the 'sigma' in the constraint term

	// Get the name of the 'sigma' for the gaussian
	// (I don't know of a way of doing RooGaussian::GetSigma() )
	// For alpha's, the sigma points to a global RooConstVar
	// with the name "1"
	// For gamma_stat_*, the sigma is named *_sigma
	std::string sigmaName = "";
	if( parameter.find("alpha_")!=std::string::npos ) {
	  sigmaName = "1";;
	}
	else if( parameter.find("gamma_stat_")!=std::string::npos ) {
	  sigmaName = parameter + "_sigma";
	}

	// Get the sigma and its value
	RooAbsReal* sigmaVar = dynamic_cast<RooAbsReal*>(constraintTerm->findServer(sigmaName.c_str()));
	if( sigmaVar==NULL ) {
	  std::cout << "Error: Failed to find the 'sigma' node: " << sigmaName
		    << " in the RooGaussian: " << constraintTerm->GetName() << std::endl;
	  throw hf_exc();
	}
	// If we find the uncertainty:
	sigma = sigmaVar->getVal();
      }
      else if( ConstraintType == "RooPoisson" ){
	// Poisson errors are given by inverting: tau = 1 / (sigma*sigma)
	std::string tauName = "nom_" + parameter;
	RooAbsReal* tauVar = dynamic_cast<RooAbsReal*>( constraintTerm->findServer(tauName.c_str()) );
	if( tauVar==NULL ) {
	  std::cout << "Error: Failed to find the nominal 'tau' node: " << tauName
		    << " for the RooPoisson: " << constraintTerm->GetName() << std::endl;
	  throw hf_exc();
	}
	double tau_val = tauVar->getVal();
	sigma = 1.0 / TMath::Sqrt( tau_val );
      } 
      else {
	std::cout << "Error: Encountered unknown constraint type for Stat Uncertainties: " 
		  << ConstraintType << std::endl;
	throw hf_exc();
      }

      return sigma;

    }

    void HistFactoryNavigation::ReplaceNode(const std::string& ToReplace, RooAbsArg* ReplaceWith) {

      // First, check that the node to replace is actually a node:
      RooAbsArg* nodeToReplace = findChild(ToReplace, fModel);
      if( nodeToReplace==NULL ) {
        std::cout << "Error: Cannot replace node: " << ToReplace
            << " because this node wasn't found in: " << fModel->GetName()
            << std::endl;
        throw hf_exc();
      }

      // Now that we have the node we want to replace, we have to 
      // get its parent node

      // Do this by looping over the clients and replacing their servers
      // (NOTE: This happens for ALL clients across the pdf)
      for (auto client : nodeToReplace->clients()) {

        // Check if this client is a member of our pdf
        // (We probably don't want to mess with clients
        // if they aren't...)
        if( findChild(client->GetName(), fModel) == nullptr) continue;

        // Now, do the replacement:
        bool valueProp=false;
        bool shapeProp=false;
        client->replaceServer( *nodeToReplace, *ReplaceWith, valueProp, shapeProp );
        std::cout << "Replaced: " << ToReplace << " with: " << ReplaceWith->GetName()
		          << " in node: " << client->GetName() << std::endl;

      }

      return;

    }


    void HistFactoryNavigation::PrintSampleComponents(const std::string& channel, 
						      const std::string& sample) {
      
      // Get the Sample Node
      RooAbsReal* sampleNode = SampleFunction(channel, sample);

      // Get the observables for this channel
      RooArgList observable_list( *GetObservableSet(channel) );

      // Make the total histogram for this sample
      std::string total_Name = sampleNode->GetName();
      TH1* total_hist= MakeHistFromRooFunction( sampleNode, observable_list, total_Name + "_tmp");
      unsigned int num_bins = total_hist->GetNbinsX()*total_hist->GetNbinsY()*total_hist->GetNbinsZ();

      RooArgSet components;
      
      // Let's see what it is...
      int label_print_width = 30;  
      int bin_print_width = 12;
      if( strcmp(sampleNode->ClassName(),"RooProduct")==0){
	RooProduct* prod = dynamic_cast<RooProduct*>(sampleNode);
	components.add( _GetAllProducts(prod) );
      }
      else {
	components.add(*sampleNode);
      }
      
      /////// NODE SIZE
      {
	TIterator* itr = components.createIterator();
	RooAbsArg* arg = NULL;
	while( (arg=(RooAbsArg*)itr->Next()) ) {
	  RooAbsReal* component = dynamic_cast<RooAbsReal*>(arg);
	  std::string NodeName = component->GetName();
	  label_print_width = TMath::Max(label_print_width, (int)NodeName.size()+2);
	}
      }

      // Now, loop over the components and print them out:
      std::cout << std::endl;
      std::cout << "Channel: " << channel << " Sample: " << sample << std::endl;
      std::cout << std::setw(label_print_width) << "Factor";

      for(unsigned int i=0; i < num_bins; ++i) {
	if( _minBinToPrint != -1 && (int)i < _minBinToPrint) continue;
	if( _maxBinToPrint != -1 && (int)i > _maxBinToPrint) break;
	std::stringstream sstr;
	sstr << "Bin" << i;
	std::cout << std::setw(bin_print_width) << sstr.str();
      }
      std::cout << std::endl;

      TIterator* itr = components.createIterator();
      RooAbsArg* arg = NULL;
      while( (arg=(RooAbsArg*)itr->Next()) ) {
	RooAbsReal* component = dynamic_cast<RooAbsReal*>(arg);
	std::string NodeName = component->GetName();

	// Make a histogram for this node	
	// Do some horrible things to prevent some really
	// annoying messages from being printed
	RooFit::MsgLevel levelBefore = RooMsgService::instance().globalKillBelow();
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	TH1* hist=NULL;
	try {
	  hist = MakeHistFromRooFunction( component, observable_list, NodeName+"_tmp");
	} catch(...) {
	  RooMsgService::instance().setGlobalKillBelow(levelBefore);
	  throw;
	}
	RooMsgService::instance().setGlobalKillBelow(levelBefore);

	// Print the hist
	std::cout << std::setw(label_print_width) << NodeName;

	// Print the Histogram
	PrintMultiDimHist(hist, bin_print_width);
	delete hist;
      }
      /////
      std::string line_break;
      int high_bin = _maxBinToPrint==-1 ? num_bins : TMath::Min(_maxBinToPrint, (int)num_bins);
      int low_bin = _minBinToPrint==-1 ? 1 : _minBinToPrint;
      int break_length = (high_bin - low_bin + 1) * bin_print_width;
      break_length += label_print_width;
      for(int i = 0; i < break_length; ++i) {
	line_break += "=";
      }
      std::cout << line_break << std::endl;

      std::cout << std::setw(label_print_width) << "TOTAL:";
      PrintMultiDimHist(total_hist, bin_print_width);
      /*
      for(unsigned int i = 0; i < num_bins; ++i) {
	if( _minBinToPrint != -1 && (int)i < _minBinToPrint) continue;
	if( _maxBinToPrint != -1 && (int)i > _maxBinToPrint) break;
	std::cout << std::setw(bin_print_width) << total_hist->GetBinContent(i+1);
      }
      std::cout << std::endl << std::endl;
      */
      delete total_hist;
      
      return;

    }


    TH1* HistFactoryNavigation::MakeHistFromRooFunction( RooAbsReal* func, RooArgList vars, 
							 std::string name ) {

      // Turn a RooAbsReal* into a TH1* based 
      // on a template histogram.  
      // The 'vars' arg list defines the x (and y and z variables)
      // Loop over the bins of the Template,
      // find the bin centers, 
      // Scan the input Var over those bin centers,
      // and use the value of the function
      // to make the new histogram

      // Make the new histogram
      // Cone and empty the template
      //      TH1* hist = (TH1*) histTemplate.Clone( name.c_str() );

      int dim = vars.getSize();

      TH1* hist=NULL;

      if( dim==1 ) {
	RooRealVar* varX = (RooRealVar*) vars.at(0);
	hist = func->createHistogram( name.c_str(),*varX, RooFit::Binning(varX->getBinning()), RooFit::Scaling(false) );
      }
      else if( dim==2 ) {
	RooRealVar* varX = (RooRealVar*) vars.at(0);
	RooRealVar* varY = (RooRealVar*) vars.at(1);
	hist = func->createHistogram( name.c_str(),*varX, RooFit::Binning(varX->getBinning()), RooFit::Scaling(false),
				      RooFit::YVar(*varY, RooFit::Binning(varY->getBinning())) );
      }
      else if( dim==3 ) {
	RooRealVar* varX = (RooRealVar*) vars.at(0);
	RooRealVar* varY = (RooRealVar*) vars.at(1);
	RooRealVar* varZ = (RooRealVar*) vars.at(2);
	hist = func->createHistogram( name.c_str(),*varX, RooFit::Binning(varX->getBinning()), RooFit::Scaling(false),
				      RooFit::YVar(*varY, RooFit::Binning(varY->getBinning())),
				      RooFit::YVar(*varZ, RooFit::Binning(varZ->getBinning())) );
      }
      else {
	std::cout << "Error: To Create Histogram from RooAbsReal function, Dimension must be 1, 2, or 3" << std::endl;
	throw hf_exc();
      }

      return hist;
    }

    // A simple wrapper to use a ModelConfig
    void HistFactoryNavigation::_GetNodes(ModelConfig* mc) {
      RooAbsPdf* modelPdf = mc->GetPdf();
      const RooArgSet* observables = mc->GetObservables();
      _GetNodes(modelPdf, observables);
    }


    void HistFactoryNavigation::SetConstant(const std::string& regExpr, bool constant) {

      // Regex FTW
      
      TString RegexTString(regExpr);
      TRegexp theRegExpr(RegexTString);

      // Now, loop over all variables and 
      // set the constant as 

      // Get the list of parameters
      RooArgSet* params = fModel->getParameters(*fObservables);
      
      std::cout << std::endl;

      // Create the title row
      std::cout << std::setw(30) << "Parameter";
      std::cout << std::setw(15) << "Value"
		<< std::setw(15) << "Error Low" 
		<< std::setw(15) << "Error High"
		<< std::endl;
      
      // Loop over the parameters and print their values, etc
      TIterator* paramItr = params->createIterator();
      RooRealVar* param = NULL;
      while( (param=(RooRealVar*)paramItr->Next()) ) {

	std::string ParamName = param->GetName();
	TString ParamNameTString(ParamName);

	// Use the Regex to skip all parameters that don't match
	//if( theRegExpr.Index(ParamNameTString, ParamName.size()) == -1 ) continue;
	Ssiz_t dummy;
	if( theRegExpr.Index(ParamNameTString, &dummy) == -1 ) continue;
	
	param->setConstant( constant );
	std::cout << "Setting param: " << ParamName << " constant" 
		  << " (matches regex: " << regExpr << ")" << std::endl;
      }
    }

    RooRealVar* HistFactoryNavigation::var(const std::string& varName) const {
      
      RooAbsArg* arg = findChild(varName, fModel);
      if( !arg ) return NULL;

      RooRealVar* var_obj = dynamic_cast<RooRealVar*>(arg);
      return var_obj;

    }

    /*
      void HistFactoryNavigation::AddChannel(const std::string& channel, RooAbsPdf* pdf, 
      RooDataSet* data=NULL) {

      }
    */

  } // namespace HistFactory
} // namespace RooStats




