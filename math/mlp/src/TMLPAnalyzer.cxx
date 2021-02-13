// @(#)root/mlp:$Id$
// Author: Christophe.Delaere@cern.ch   25/04/04

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TMLPAnalyzer

This utility class contains a set of tests usefull when developing
a neural network.
It allows you to check for unneeded variables, and to control
the network structure.

*/

#include "TROOT.h"
#include "TSynapse.h"
#include "TNeuron.h"
#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TEventList.h"
#include "TH1D.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TVirtualPad.h"
#include "TRegexp.h"
#include "TMath.h"
#include "snprintf.h"
#include <iostream>
#include <cstdlib>

ClassImp(TMLPAnalyzer);

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TMLPAnalyzer::~TMLPAnalyzer()
{
   delete fAnalysisTree;
   delete fIOTree;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the number of layers.

Int_t TMLPAnalyzer::GetLayers()
{
   TString fStructure = fNetwork->GetStructure();
   return fStructure.CountChar(':')+1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the number of neurons in given layer.

Int_t TMLPAnalyzer::GetNeurons(Int_t layer)
{
   if(layer==1) {
      TString fStructure = fNetwork->GetStructure();
      TString input      = TString(fStructure(0, fStructure.First(':')));
      return input.CountChar(',')+1;
   }
   else if(layer==GetLayers()) {
      TString fStructure = fNetwork->GetStructure();
      TString output = TString(fStructure(fStructure.Last(':') + 1,
                               fStructure.Length() - fStructure.Last(':')));
      return output.CountChar(',')+1;
   }
   else {
      Int_t cnt=1;
      TString fStructure = fNetwork->GetStructure();
      TString hidden = TString(fStructure(fStructure.First(':') + 1,
                               fStructure.Last(':') - fStructure.First(':') - 1));
      Int_t beg = 0;
      Int_t end = hidden.Index(":", beg + 1);
      Int_t num = 0;
      while (end != -1) {
         num = atoi(TString(hidden(beg, end - beg)).Data());
         cnt++;
         beg = end + 1;
         end = hidden.Index(":", beg + 1);
         if(layer==cnt) return num;
      }
      num = atoi(TString(hidden(beg, hidden.Length() - beg)).Data());
      cnt++;
      if(layer==cnt) return num;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the formula used as input for neuron (idx) in
/// the first layer.

TString TMLPAnalyzer::GetNeuronFormula(Int_t idx)
{
   TString fStructure = fNetwork->GetStructure();
   TString input      = TString(fStructure(0, fStructure.First(':')));
   Int_t beg = 0;
   Int_t end = input.Index(",", beg + 1);
   TString brName;
   Int_t cnt = 0;
   while (end != -1) {
      brName = TString(input(beg, end - beg));
      if (brName[0]=='@')
         brName = brName(1,brName.Length()-1);
      beg = end + 1;
      end = input.Index(",", beg + 1);
      if(cnt==idx) return brName;
      cnt++;
   }
   brName = TString(input(beg, input.Length() - beg));
   if (brName[0]=='@')
      brName = brName(1,brName.Length()-1);
   return brName;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the name of any neuron from the input layer

const char* TMLPAnalyzer::GetInputNeuronTitle(Int_t in)
{
   TNeuron* neuron=(TNeuron*)fNetwork->fFirstLayer[in];
   return neuron ? neuron->GetName() : "NO SUCH NEURON";
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the name of any neuron from the output layer

const char* TMLPAnalyzer::GetOutputNeuronTitle(Int_t out)
{
   TNeuron* neuron=(TNeuron*)fNetwork->fLastLayer[out];
   return neuron ? neuron->GetName() : "NO SUCH NEURON";
}

////////////////////////////////////////////////////////////////////////////////
/// Gives some information about the network in the terminal.

void TMLPAnalyzer::CheckNetwork()
{
   TString fStructure = fNetwork->GetStructure();
   std::cout << "Network with structure: " << fStructure.Data() << std::endl;
   std::cout << "inputs with low values in the differences plot may not be needed" << std::endl;
   // Checks if some input variable is not needed
   char var[64], sel[64];
   for (Int_t i = 0; i < GetNeurons(1); i++) {
      snprintf(var,64,"diff>>tmp%d",i);
      snprintf(sel,64,"inNeuron==%d",i);
      fAnalysisTree->Draw(var, sel, "goff");
      TH1F* tmp = (TH1F*)gDirectory->Get(Form("tmp%d",i));
      if (!tmp) continue;
      std::cout << GetInputNeuronTitle(i)
           << " -> " << tmp->GetMean()
           << " +/- " << tmp->GetRMS() << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Collect information about what is useful in the network.
/// This method has to be called first when analyzing a network.
/// Fills the two analysis trees.

void TMLPAnalyzer::GatherInformations()
{
   Double_t shift = 0.1;
   TTree* data = fNetwork->fData;
   TEventList* test = fNetwork->fTest;
   Int_t nEvents = test->GetN();
   Int_t nn = GetNeurons(1);
   Double_t* params = new Double_t[nn];
   Double_t* rms    = new Double_t[nn];
   TTreeFormula** formulas = new TTreeFormula*[nn];
   Int_t* index = new Int_t[nn];
   TString formula;
   TRegexp re("{[0-9]+}$");
   Ssiz_t len = formula.Length();
   Ssiz_t pos = -1;
   Int_t i(0), j(0), k(0), l(0);
   for(i=0; i<nn; i++){
      formula = GetNeuronFormula(i);
      pos = re.Index(formula,&len);
      if(pos==-1 || len<3) {
         formulas[i] = new TTreeFormula(Form("NF%lu",(ULong_t)this),formula,data);
         index[i] = 0;
      }
      else {
         TString newformula(formula,pos);
         TString val = formula(pos+1,len-2);
         formulas[i] = new TTreeFormula(Form("NF%lu",(ULong_t)this),newformula,data);
         formula = newformula;
         index[i] = val.Atoi();
      }
      TH1D tmp("tmpb", "tmpb", 1, -FLT_MAX, FLT_MAX);
      data->Draw(Form("%s>>tmpb",formula.Data()),"","goff");
      rms[i]  = tmp.GetRMS();
   }
   Int_t inNeuron = 0;
   Double_t diff = 0.;
   if(fAnalysisTree) delete fAnalysisTree;
   fAnalysisTree = new TTree("result","analysis");
   fAnalysisTree->SetDirectory(0);
   fAnalysisTree->Branch("inNeuron",&inNeuron,"inNeuron/I");
   fAnalysisTree->Branch("diff",&diff,"diff/D");
   Int_t numOutNodes=GetNeurons(GetLayers());
   Double_t *outVal=new Double_t[numOutNodes];
   Double_t *trueVal=new Double_t[numOutNodes];

   delete fIOTree;
   fIOTree=new TTree("MLP_iotree","MLP_iotree");
   fIOTree->SetDirectory(0);
   TString leaflist;
   for (i=0; i<nn; i++)
      leaflist+=Form("In%d/D:",i);
   leaflist.Remove(leaflist.Length()-1);
   fIOTree->Branch("In", params, leaflist);

   leaflist="";
   for (i=0; i<numOutNodes; i++)
      leaflist+=Form("Out%d/D:",i);
   leaflist.Remove(leaflist.Length()-1);
   fIOTree->Branch("Out", outVal, leaflist);

   leaflist="";
   for (i=0; i<numOutNodes; i++)
      leaflist+=Form("True%d/D:",i);
   leaflist.Remove(leaflist.Length()-1);
   fIOTree->Branch("True", trueVal, leaflist);
   Double_t v1 = 0.;
   Double_t v2 = 0.;
   // Loop on the events in the test sample
   for(j=0; j< nEvents; j++) {
      fNetwork->GetEntry(test->GetEntry(j));
      // Loop on the neurons to evaluate
      for(k=0; k<GetNeurons(1); k++) {
         params[k] = formulas[k]->EvalInstance(index[k]);
      }
      for(k=0; k<GetNeurons(GetLayers()); k++) {
         outVal[k] = fNetwork->Evaluate(k,params);
         trueVal[k] = ((TNeuron*)fNetwork->fLastLayer[k])->GetBranch();
      }
      fIOTree->Fill();

      // Loop on the input neurons
      for (i = 0; i < GetNeurons(1); i++) {
         inNeuron = i;
         diff = 0;
         // Loop on the neurons in the output layer
         for(l=0; l<GetNeurons(GetLayers()); l++){
            params[i] += shift*rms[i];
            v1 = fNetwork->Evaluate(l,params);
            params[i] -= 2*shift*rms[i];
            v2 = fNetwork->Evaluate(l,params);
            diff += (v1-v2)*(v1-v2);
            // reset to original value
            params[i] += shift*rms[i];
         }
         diff = TMath::Sqrt(diff);
         fAnalysisTree->Fill();
      }
   }
   delete[] params;
   delete[] rms;
   delete[] outVal;
   delete[] trueVal;
   delete[] index;
   for(i=0; i<GetNeurons(1); i++) delete formulas[i];
   delete [] formulas;
   fAnalysisTree->ResetBranchAddresses();
   fIOTree->ResetBranchAddresses();
}

////////////////////////////////////////////////////////////////////////////////
/// Draws the distribution (on the test sample) of the
/// impact on the network output of a small variation of
/// the ith input.

void TMLPAnalyzer::DrawDInput(Int_t i)
{
   char sel[64];
   snprintf(sel,64, "inNeuron==%d", i);
   fAnalysisTree->Draw("diff", sel);
}

////////////////////////////////////////////////////////////////////////////////
/// Draws the distribution (on the test sample) of the
/// impact on the network output of a small variation of
/// each input.
/// DrawDInputs() draws something that approximates the distribution of the
/// derivative of the NN w.r.t. each input. That quantity is recognized as
/// one of the measures to determine key quantities in the network.
///
/// What is done is to vary one input around its nominal value and to see
/// how the NN changes. This is done for each entry in the sample and produces
/// a distribution.
///
/// What you can learn from that is:
/// - is variable a really useful, or is my network insensitive to it ?
/// - is there any risk of big systematic ? Is the network extremely sensitive
///   to small variations of any of my inputs ?
///
/// As you might understand, this is to be considered with care and can serve
/// as input for an "educated guess" when optimizing the network.

void TMLPAnalyzer::DrawDInputs()
{
   THStack* stack  = new THStack("differences","differences (impact of variables on ANN)");
   TLegend* legend = new TLegend(0.75,0.75,0.95,0.95);
   TH1F* tmp = 0;
   char var[64], sel[64];
   for(Int_t i = 0; i < GetNeurons(1); i++) {
      snprintf(var,64, "diff>>tmp%d", i);
      snprintf(sel,64, "inNeuron==%d", i);
      fAnalysisTree->Draw(var, sel, "goff");
      tmp = (TH1F*)gDirectory->Get(Form("tmp%d",i));
      tmp->SetDirectory(0);
      tmp->SetLineColor(i+1);
      stack->Add(tmp);
      legend->AddEntry(tmp,GetInputNeuronTitle(i),"l");
   }
   stack->Draw("nostack");
   legend->Draw();
   gPad->SetLogy();
}

////////////////////////////////////////////////////////////////////////////////
/// Draws the distribution of the neural network (using ith neuron).
/// Two distributions are drawn, for events passing respectively the "signal"
/// and "background" cuts. Only the test sample is used.

void TMLPAnalyzer::DrawNetwork(Int_t neuron, const char* signal, const char* bg)
{
   TTree* data = fNetwork->fData;
   TEventList* test = fNetwork->fTest;
   TEventList* current = data->GetEventList();
   data->SetEventList(test);
   THStack* stack = new THStack("__NNout_TMLPA",Form("Neural net output (neuron %d)",neuron));
   TH1F *bgh  = new TH1F("__bgh_TMLPA", "NN output", 50, -0.5, 1.5);
   TH1F *sigh = new TH1F("__sigh_TMLPA", "NN output", 50, -0.5, 1.5);
   bgh->SetDirectory(0);
   sigh->SetDirectory(0);
   Int_t nEvents = 0;
   Int_t j=0;
   // build event lists for signal and background
   TEventList* signal_list = new TEventList("__tmpSig_MLPA");
   TEventList* bg_list     = new TEventList("__tmpBkg_MLPA");
   data->Draw(">>__tmpSig_MLPA",signal,"goff");
   data->Draw(">>__tmpBkg_MLPA",bg,"goff");

   // fill the background
   nEvents = bg_list->GetN();
   for(j=0; j< nEvents; j++) {
      bgh->Fill(fNetwork->Result(bg_list->GetEntry(j),neuron));
   }
   // fill the signal
   nEvents = signal_list->GetN();
   for(j=0; j< nEvents; j++) {
      sigh->Fill(fNetwork->Result(signal_list->GetEntry(j),neuron));
   }
   // draws the result
   bgh->SetLineColor(kBlue);
   bgh->SetFillStyle(3008);
   bgh->SetFillColor(kBlue);
   sigh->SetLineColor(kRed);
   sigh->SetFillStyle(3003);
   sigh->SetFillColor(kRed);
   bgh->SetStats(0);
   sigh->SetStats(0);
   stack->Add(bgh);
   stack->Add(sigh);
   TLegend *legend = new TLegend(.75, .80, .95, .95);
   legend->AddEntry(bgh, "Background");
   legend->AddEntry(sigh,"Signal");
   stack->Draw("nostack");
   legend->Draw();
   // restore the default event list
   data->SetEventList(current);
   delete signal_list;
   delete bg_list;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a profile of the difference of the MLP output minus the
/// true value for a given output node outnode, vs the true value for
/// outnode, for all test data events. This method is mainly useful
/// when doing regression analysis with the MLP (i.e. not classification,
/// but continuous truth values).
/// The resulting TProfile histogram is returned.
/// It is not drawn if option "goff" is specified.
/// Options are passed to TProfile::Draw

TProfile* TMLPAnalyzer::DrawTruthDeviation(Int_t outnode /*=0*/,
                                           Option_t *option /*=""*/)
{
   if (!fIOTree) GatherInformations();
   TString pipehist=Form("MLP_truthdev_%d",outnode);
   TString drawline;
   drawline.Form("Out.Out%d-True.True%d:True.True%d>>",
                 outnode, outnode, outnode);
   fIOTree->Draw(drawline+pipehist+"(20)", "", "goff prof");
   TProfile* h=(TProfile*)gDirectory->Get(pipehist);
   h->SetDirectory(0);
   const char* title=GetOutputNeuronTitle(outnode);
   if (title) {
      h->SetTitle(Form("#Delta(output - truth) vs. truth for %s",
                      title));
      h->GetXaxis()->SetTitle(title);
      h->GetYaxis()->SetTitle(Form("#Delta(output - truth) for %s", title));
   }
   if (!strstr(option,"goff"))
      h->Draw();
   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates TProfiles of the difference of the MLP output minus the
/// true value vs the true value, one for each output, filled with the
/// test data events. This method is mainly useful when doing regression
/// analysis with the MLP (i.e. not classification, but continuous truth
/// values).
/// The returned THStack contains all the TProfiles. It is drawn unless
/// the option "goff" is specified.
/// Options are passed to TProfile::Draw.

THStack* TMLPAnalyzer::DrawTruthDeviations(Option_t *option /*=""*/)
{
   THStack *hs=new THStack("MLP_TruthDeviation",
                           "Deviation of MLP output from truth");

   // leg!=0 means we're drawing
   TLegend *leg=0;
   if (!option || !strstr(option,"goff"))
      leg=new TLegend(.4,.85,.95,.95,"#Delta(output - truth) vs. truth for:");

   const char* xAxisTitle=0;

   // create profile for each input neuron,
   // adding them into the THStack and the TLegend
   for (Int_t outnode=0; outnode<GetNeurons(GetLayers()); outnode++) {
      TProfile* h=DrawTruthDeviation(outnode, "goff");
      h->SetLineColor(1+outnode);
      hs->Add(h, option);
      if (leg) leg->AddEntry(h,GetOutputNeuronTitle(outnode));
      if (!outnode)
         // Xaxis title is the same for all, extract it from the first one.
         xAxisTitle=h->GetXaxis()->GetTitle();
   }

   if (leg) {
      hs->Draw("nostack");
      leg->Draw();
      // gotta draw before accessing the axes
      hs->GetXaxis()->SetTitle(xAxisTitle);
      hs->GetYaxis()->SetTitle("#Delta(output - truth)");
   }

   return hs;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a profile of the difference of the MLP output outnode minus
/// the true value of outnode vs the input value innode, for all test
/// data events.
/// The resulting TProfile histogram is returned.
/// It is not drawn if option "goff" is specified.
/// Options are passed to TProfile::Draw

TProfile* TMLPAnalyzer::DrawTruthDeviationInOut(Int_t innode,
                                                Int_t outnode /*=0*/,
                                                Option_t *option /*=""*/)
{
   if (!fIOTree) GatherInformations();
   TString pipehist=Form("MLP_truthdev_i%d_o%d", innode, outnode);
   TString drawline;
   drawline.Form("Out.Out%d-True.True%d:In.In%d>>",
                 outnode, outnode, innode);
   fIOTree->Draw(drawline+pipehist+"(50)", "", "goff prof");
   TProfile* h=(TProfile*)gROOT->FindObject(pipehist);
   h->SetDirectory(0);
   const char* titleInNeuron=GetInputNeuronTitle(innode);
   const char* titleOutNeuron=GetOutputNeuronTitle(outnode);
   h->SetTitle(Form("#Delta(output - truth) of %s vs. input %s",
                    titleOutNeuron, titleInNeuron));
   h->GetXaxis()->SetTitle(Form("%s", titleInNeuron));
   h->GetYaxis()->SetTitle(Form("#Delta(output - truth) for %s",
                                titleOutNeuron));
   if (!strstr(option,"goff"))
      h->Draw(option);
   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a profile of the difference of the MLP output outnode minus the
/// true value of outnode vs the input value, stacked for all inputs, for
/// all test data events.
/// The returned THStack contains all the TProfiles. It is drawn unless
/// the option "goff" is specified.
/// Options are passed to TProfile::Draw.

THStack* TMLPAnalyzer::DrawTruthDeviationInsOut(Int_t outnode /*=0*/,
                                                Option_t *option /*=""*/)
{
   TString sName;
   sName.Form("MLP_TruthDeviationIO_%d", outnode);
   const char* outputNodeTitle=GetOutputNeuronTitle(outnode);
   THStack *hs=new THStack(sName,
                           Form("Deviation of MLP output %s from truth",
                                outputNodeTitle));

   // leg!=0 means we're drawing.
   TLegend *leg=0;
   if (!option || !strstr(option,"goff"))
      leg=new TLegend(.4,.75,.95,.95,
                      Form("#Delta(output - truth) of %s vs. input for:",
                           outputNodeTitle));

   // create profile for each input neuron,
   // adding them into the THStack and the TLegend
   Int_t numInNodes=GetNeurons(1);
   Int_t innode=0;
   for (innode=0; innode<numInNodes; innode++) {
      TProfile* h=DrawTruthDeviationInOut(innode, outnode, "goff");
      h->SetLineColor(1+innode);
      hs->Add(h, option);
      if (leg) leg->AddEntry(h,h->GetXaxis()->GetTitle());
   }

   if (leg) {
      hs->Draw("nostack");
      leg->Draw();
      // gotta draw before accessing the axes
      hs->GetXaxis()->SetTitle("Input value");
      hs->GetYaxis()->SetTitle(Form("#Delta(output - truth) for %s",
                                 outputNodeTitle));
   }

   return hs;
}
