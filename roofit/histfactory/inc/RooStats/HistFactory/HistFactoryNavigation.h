#ifndef INCLUDE_HISTFACTORYNAVIGATION_H
#define INCLUDE_HISTFACTORYNAVIGATION_H

#include <map>
#include <vector>
#include <string>

#include "TH1.h"
#include "THStack.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooProduct.h"
#include "RooStats/HistFactory/Measurement.h"

namespace RooStats {
  class ModelConfig;
  namespace HistFactory {

    class HistFactoryNavigation {

    public:

      /// Initialze based on an already-created HistFactory Model
      HistFactoryNavigation(ModelConfig* mc);
      HistFactoryNavigation(const std::string& File, 
			    const std::string& WorkspaceName="combined",
			    const std::string& ModelConfigName="ModelConfig");
      HistFactoryNavigation(RooAbsPdf* model, RooArgSet* observables);

      virtual ~HistFactoryNavigation() {} 

      /// Should pretty print all channels and the current values 
      void PrintState();
      /// Should pretty print this and the current values
      void PrintState(const std::string& channel);

      /// Print the current values and errors of pdf parameters
      void PrintParameters(bool IncludeConstantParams=false);

      /// Print parameters that effect a particular channel
      void PrintChannelParameters(const std::string& channel, 
				  bool IncludeConstantParams=false);

      /// Print parameters that effect a particular sample
      void PrintSampleParameters(const std::string& channel, const std::string& sample, 
				 bool IncludeConstantParams=false);

      /// Print the different components that make up a sample
      /// (NormFactors, Statistical Uncertainties, Interpolation, etc)
      void PrintSampleComponents(const std::string& channel, const std::string& sample);

      /// Print a "HistFactory style" RooDataSet in a readable way
      void PrintDataSet(RooDataSet* data, const std::string& channel="");

      /// Print the model and the data, comparing channel by channel
      void PrintModelAndData(RooDataSet* data);

      /// The value of the ith bin for the total in that channel
      double GetBinValue(int bin, const std::string& channel);  
      /// The value of the ith bin for that sample and channel 
      double GetBinValue(int bin, const std::string& channel, const std::string& sample);  

      /// The (current) histogram for that sample 
      /// This includes all parameters and interpolation
      TH1* GetSampleHist(const std::string& channel, 
			 const std::string& sample, const std::string& name="");  

      /// Get the total channel histogram for this channel
      TH1* GetChannelHist(const std::string& channel, const std::string& name="");

      /// Get a histogram from the dataset for this channel
      TH1* GetDataHist(RooDataSet* data, const std::string& channel, const std::string& name="");

      /// Get a stack of all samples in a channel
      THStack* GetChannelStack(const std::string& channel, const std::string& name="");

      /// Draw a stack of the channel, and include data if the pointer is supplied
      void DrawChannel(const std::string& channel, RooDataSet* data=NULL);

      /// Get the RooAbsReal function for a given sample in a given channel
      RooAbsReal* SampleFunction(const std::string& channel, const std::string& sample);

      /// Get the set of observables for a given channel
      RooArgSet* GetObservableSet(const std::string& channel);

      /// Get the constraint term for a given systematic (alpha or gamma)
      RooAbsReal* GetConstraintTerm(const std::string& parameter);

      /// Get the uncertainty based on the constraint term for a given systematic
      double GetConstraintUncertainty(const std::string& parameter);

      /// Find a node in the pdf and replace it with a new node
      /// These nodes can be functions, pdf's, RooRealVar's, etc
      /// Will do minimial checking to make sure the replacement makes sense
      void ReplaceNode(const std::string& ToReplace, RooAbsArg* ReplaceWith);

      // Set any RooRealVar's const (or not const) if they match
      // the supplied regular expression
      void SetConstant(const std::string& regExpr=".*", bool constant=true);

      void SetMaxBinToPrint(int max) { _maxBinToPrint = max; }
      int GetMaxBinToPrint() const { return _maxBinToPrint; }

      void SetMinBinToPrint(int min) { _minBinToPrint = min; }
      int GetMinBinToPrint() const { return _minBinToPrint; }

      /// Get the model for this channel
      RooAbsPdf* GetModel() const { return fModel; }

      //
      RooAbsPdf* GetChannelPdf(const std::string& channel);


      std::vector< std::string > GetChannelSampleList(const std::string& channel);

      // Return the RooRealVar by the same name used in the model
      // If not found, return NULL
      RooRealVar* var(const std::string& varName) const;
      
      /*
      // Add a channel to the pdf
      // Combine the data if it is provided
      void AddChannel(const std::string& channel, RooAbsPdf* pdf, RooDataSet* data=NULL);
      */

      /*
      // Add a constraint term to the pdf
      // This method requires that the pdf NOT be simultaneous
      void AddConstraintTerm(RooAbsArg* constraintTerm);

      // Add a constraint term to the pdf of a particular channel
      // This method requires that the pdf be simultaneous 
      // OR that the channel string match the channel that the pdf represents
      void AddConstraintTerm(RooAbsArg* constraintTerm, const std::string& channel);
      */

    protected:

      /// Set the title and bin widths
      void SetPrintWidths(const std::string& channel);

      /// Fetch the node information for the pdf in question, and
      /// save it in the varous collections in this class
      void _GetNodes(ModelConfig* mc);
      void _GetNodes(RooAbsPdf* model, const RooArgSet* observables);

      /// Print a histogram's contents to the screen
      /// void PrettyPrintHistogram(TH1* hist);
      void PrintMultiDimHist(TH1* hist, int bin_print_width);

      /// Make a histogram from a funciton
      /// Edit so it can take a RooArgSet of parameters
      TH1* MakeHistFromRooFunction( RooAbsReal* func, RooArgList vars, std::string name="Hist" );

      /// Get a map of sample names to their functions for a particular channel
      std::map< std::string, RooAbsReal*> GetSampleFunctionMap(const std::string& channel);

    private:
      
      /// The HistFactory Pdf Pointer
      RooAbsPdf* fModel;

      /// The observables
      RooArgSet* fObservables;

      int _minBinToPrint;
      int _maxBinToPrint;

      int _label_print_width;
      int _bin_print_width;

      /// The list of channels
      std::vector<std::string> fChannelNameVec;

      /// Map of channel names to their full pdf's
      std::map< std::string, RooAbsPdf* > fChannelPdfMap;  

      /// Map of channel names to pdf without constraint
      std::map< std::string, RooAbsPdf* > fChannelSumNodeMap;  
      
      /// Map of channel names to their set of ovservables
      std::map< std::string, RooArgSet*> fChannelObservMap;

      /// Map of Map of Channel, Sample names to Function Nodes
      /// Used by doing: fChannelSampleFunctionMap["MyChannel"]["MySample"]
      std::map< std::string, std::map< std::string, RooAbsReal*> > fChannelSampleFunctionMap;

      /// Internal method implementation of finding a daughter node
      /// from a parent node (looping over all generations)
      RooAbsArg* findChild(const std::string& name, RooAbsReal* parent) const;

      /// Recursively get all products of products
      RooArgSet _GetAllProducts(RooProduct* node);
      

    protected:
      ClassDef(RooStats::HistFactory::HistFactoryNavigation,2)

    };

  }
}
#endif // INCLUDE_HISTFACTORYNAVIGATION_H
