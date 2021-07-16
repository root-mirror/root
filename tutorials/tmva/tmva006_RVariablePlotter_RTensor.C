//
/// \file
/// \ingroup tutorial_tmva
/// Authors: Simone Azeglio, Lorenzo Moneta, Sitong An, Stefan Wunsch
/// This tutorial shows how to plot the distribution of variables from RTensor Data
///


#include "TMVA/RVariablePlotter.h"
#include "TMVA/tmvaglob.h"

using namespace TMVA::Experimental;


void tmva006_RVariablePlotter_Tensor()
{
    // Initialize ROOT dataframes from signal and background datasets
    const std::string filename = "http://root.cern.ch/files/tmva_class_example.root";
    ROOT::RDataFrame sig1("TreeS", filename);
    ROOT::RDataFrame bkg1("TreeB", filename);
    
    // Apply transformations on the datasets to be included in the study
    auto transform_ = [](ROOT::RDF::RNode df) { return df.Define("var5", "var1 * var2"); };
    auto sig2 = transform_(sig1);
    auto bkg2 = transform_(bkg1);


    // same data as in tmva005.... , we can transform them in tensors and compare the plots
    auto sig2Tensor = AsTensor<float>(sig2);
    auto bkg2Tensor = AsTensor<float>(bkg2);
    
    // Place plots on the pads of the canvas - new columns are always added first
    //we should remember it while converting an RTensor to an RNode
    const std::vector<std::string> vars = sig2.GetColumnNames(); //{"var5", "var2", "var3", "var4", "var1"};
    
    // Create a variable plotter object giving the Tensors and the class labels.
    TMVA::RVariablePlotter plotter({sig2Tensor, bkg2Tensor}, {"Signal", "Background"});
    
    
    TCanvas *c = new TCanvas("c", "c", 1400, 800);
    c->Divide(3, 2);

    // legend vertices
    float minX = 0.7;
    float minY = 0.8;
    float maxX = 0.9;
    float maxY = 0.9;

    for (unsigned int i = 0; i < vars.size(); i++) {
        c->cd(i + 1);
        gPad->SetMargin(0.2, 0.9, 0.1, 0.9);
        c->Update();
        //gPad->SetGrid(1,1); // plotting a background grid
        plotter.DrawTensor(vars[i], vars, true);
        plotter.DrawLegend(minX, minY, maxX, maxY);
    }

}



