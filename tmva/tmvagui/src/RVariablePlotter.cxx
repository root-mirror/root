// @(#)root/tmva $Id$
// Authors: Simone Azeglio, Sitong An, Lorenzo Moneta , Stefan Wunsch

/*************************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis          *
 * Package: TMVA                                                                     *
 * Class  : RVariablePlotter                                                         *
 * Web    : http://tmva.sourceforge.net                                              *
 *                                                                                   *
 * Description:                                                                      *
 *      Variable Plotter                                                             *
 *                                                                                   *
 * Authors (alphabetical):                                                           *
 *      Simone Azeglio, University of Turin (Master Student), CERN (Summer Student)  *
 *      Lorenzo Moneta, CERN
 *      Sitong An, CERN / CMU                                                        *
 *      Stefan Wunsch                                                                *
 *                                                                                   *
 * Copyright (c) 2021:                                                               *
 *                                                                                   *
 * Redistribution and use in source and binary forms, with or without                *
 * modification, are permitted according to the terms listed in LICENSE              *
 * (http://tmva.sourceforge.net/LICENSE)                                             *
 **********************************************************************************/


/*! \class TMVA::RVariablePlotter
\ingroup TMVA
Plotting a single variable
*/

#include "TMVA/RVariablePlotter.h"
#include "TMVA/tmvaglob.h"
#include "TMVA/RTensor.hxx"


////////////////////////////////////////////////////////////////////////////////
/// constructor for RDataframe with nodes (samples) and labels

TMVA::RVariablePlotter::RVariablePlotter( const std::vector<ROOT::RDF::RNode>& nodes, const std::vector<std::string>& labels)
    : fNodes(nodes),
    fLabels(labels){
        
        if (fNodes.size() != fLabels.size())
            std::runtime_error("Number of given RDataFrame nodes does not match number of given class labels.");

        if (fNodes.size() == 0)
            std::runtime_error("Number of given RDataFrame nodes and number of given class labels cannot be zero.");
}


////////////////////////////////////////////////////////////////////////////////
/// constructor (overloading) for RTensor nodes (samples) and labels

TMVA::RVariablePlotter::RVariablePlotter( const std::vector<TMVA::Experimental::RTensor<Float_t>>& tensors, const std::vector<std::string>& labels)
    : fTensors(tensors),
    fLabels(labels)
    {
        
        if (fTensors.size() != fLabels.size())
            std::runtime_error("Number of given RTensor components does not match number of given class labels.");

        if (fTensors.size() == 0)
            std::runtime_error("Number of given RTensor components and number of given class labels cannot be zero.");
}

////////////////////////////////////////////////////////////////////////////////
/// Set style and keep existing canvas
void TMVA::RVariablePlotter::InitializeStyle(bool useTMVAStyle){
   
   // set style
   if (!useTMVAStyle) {
      gROOT->SetStyle("Plain");
      gStyle->SetOptStat(0);
      return;
   }

   TMVA::TMVAGlob::SetTMVAStyle();
}

////////////////////////////////////////////////////////////////
/// Drawing variables' plot RDataframe
void TMVA::RVariablePlotter::Draw(const std::string& variable, bool useTMVAStyle) {
   // Make histograms with TH1D
    
    TMVA::RVariablePlotter::InitializeStyle(useTMVAStyle);
    
    const auto size = fNodes.size();
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos;
    
    for (std::size_t i = 0; i < size; i++) {
        // Trigger event loop with computing the histogram
        ROOT::RDF::RResultPtr<TH1D> h = fNodes[i].Histo1D(variable);
        histos.push_back(h);
    }

   // Modify style and draw histograms
    THStack stack_;
    
    for (unsigned int i = 0; i < histos.size(); i++) {
        histos[i]->SetLineColor(i + 1);
        histos[i]->SetTitle(variable.c_str());
        histos[i]->SetStats(true);
        stack_.Add(histos[i].GetPtr());
   }
    
    auto clone = (THStack*) stack_.DrawClone("nostack");
    
    clone->SetTitle(variable.c_str());
    clone->GetXaxis()->SetTitle(variable.c_str());
    clone->GetYaxis()->SetTitle("Count");
    clone->GetYaxis()->SetTitleOffset( 1.50 );
        
}

////////////////////////////////////////////////////////////////
/// Drawing variables' plot from RTensor input
void TMVA::RVariablePlotter::DrawTensor(const std::string& variable, const std::vector<std::string>& variables, bool useTMVAStyle) {
   
    // vector of RNodes
    auto RNodeVec = TMVA::RVariablePlotter::TensorsToNodes(variables);
    
    TMVA::RVariablePlotter::InitializeStyle(useTMVAStyle);
    
    // Make histograms with TH1D
    const auto size = RNodeVec.size();
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos;
    
    for (std::size_t i = 0; i < size; i++) {
        // Trigger event loop with computing the histogram
        ROOT::RDF::RResultPtr<TH1D> h = RNodeVec[i].Histo1D(variable);
        histos.push_back(h);
        h = ROOT::RDF::RResultPtr<TH1D>();
    }

   // Modify style and draw histograms
    THStack stack;
    
    for (unsigned int i = 0; i < histos.size(); i++) {
        histos[i]->SetLineColor(i + 1);
        histos[i]->SetTitle(variable.c_str());
        histos[i]->SetStats(true);
        stack.Add(histos[i].GetPtr());
   }
    
    auto clone = (THStack*) stack.DrawClone("nostack");
    
    clone->SetTitle(variable.c_str());
    clone->GetXaxis()->SetTitle(variable.c_str());
    clone->GetYaxis()->SetTitle("Count");
    clone->GetYaxis()->SetTitleOffset( 1.50 );
        
}


////////////////////////////////////////////////////////////////////////////////
/// Drawing Legend
void TMVA::RVariablePlotter::DrawLegend(float minX = 0.8, float minY = 0.8, float maxX = 0.9, float maxY = 0.9) {
    // make Legend from TLegend
    TLegend l(minX, minY, maxX, maxY);
    std::vector<TH1D> histos(fLabels.size());

    for (unsigned int i = 0; i < fLabels.size(); i++) {
        histos[i].SetLineColor(i + 1);
        l.AddEntry(&histos[i], fLabels[i].c_str(), "l");
    }
    l.SetBorderSize(1);
    l.DrawClone();
}

////////////////////////////////////////////////////////////////////////////////
/// RTensor to RNode Converter (Static Function - Private)

// add variables argument for custom variables
ROOT::RDF::RNode TMVA::RVariablePlotter::TensorToNode(const TMVA::Experimental::RTensor<Float_t>& tensor,  const std::vector<std::string>& variables){
        
    // shape check
    if (tensor.GetShape().size() != 2)
        std::runtime_error("Number of given RTensor dimensions does not match number of given class labels.");
    
    
    // Pivot vector for RTensor to RDataFrame conversion
    //std::vector<float> vecSig;
    //ROOT::VecOps::RVec<Float_t> vecSig;
    
    // Create a RDataFrame by passing RTensor's data with a pivot vector
    //for (unsigned int i = 0; i< tensor.GetShape()[0]; i++){
      //  vecSig.push_back(tensor(i,0));}
    
    auto dfSig = ROOT::RDataFrame(tensor.GetShape()[0]).DefineSlotEntry(variables[0], [=] (unsigned, ULong64_t entry) { return tensor(entry, 0); } );
    
    std::size_t nvar = variables.size();
    // in case RTensors have a different number of columns
    if(tensor.GetShape()[1]!= variables.size())
        nvar = std::min(tensor.GetShape()[1], variables.size());
    
    for (std::size_t j = 1; j < nvar; j++){
        //vecSig = std::vector<float>();
        //vecSig = ROOT::VecOps::RVec<Float_t>();
        
        //for (unsigned int i = j*tensor.GetShape()[0]; i < (j+1)*tensor.GetShape()[0]; i++){
          //      vecSig.push_back(tensor(i,j));
            //    }
            
        dfSig = dfSig.DefineSlotEntry(variables[j],  [=] (unsigned, ULong64_t entry) { return tensor(entry,j);});
            }
    
    // RDataFrame to RNode
    auto DFNode = ROOT::RDF::RNode(dfSig);
    
    return DFNode;
}

////////////////////////////////////////////////////////////////////////////////
/// RTensor vector to RNode vector converter (Private)
std::vector<ROOT::RDF::RNode> TMVA::RVariablePlotter::TensorsToNodes(const std::vector<std::string>& variables){
    
    const auto size = fTensors.size();
    std::vector<ROOT::RDF::RNode> RNodeVec;

    // loop through tensors
    for (std::size_t k = 0; k < size; k++){
        
        RNodeVec.push_back(TensorToNode(fTensors[k], variables));
    }
    
    return RNodeVec;
    
}


