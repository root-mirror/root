// @(#)root/roostats:$Id: ProfileInspector.cxx 34109 2010-06-24 15:00:16Z moneta $

/*************************************************************************
 * Project: RooStats                                                     *
 * Package: RooFit/RooStats                                              *
 * Authors:                                                              *
 *   Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke       *
 *   Akira Shibata
 *************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


/** \class RooStats::ProfileInspector
    \ingroup Roostats

Utility class to plot conditional MLE of nuisance parameters vs. Parameters of Interest

*/


#include "RooStats/ProfileInspector.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"
#include "RooCurve.h"
#include "TAxis.h"

ClassImp(RooStats::ProfileInspector);

using namespace RooStats;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

ProfileInspector::ProfileInspector()
{
}

////////////////////////////////////////////////////////////////////////////////
/// ProfileInspector destructor

ProfileInspector::~ProfileInspector()
{
}

////////////////////////////////////////////////////////////////////////////////
///
/// This tool makes a plot of the conditional maximum likelihood estimate of the nuisance parameter
///   vs the parameter of interest
///
/// This enables you to discover if any of the nuisance parameters are behaving strangely
/// curve is the optional parameters, when used you can specify the points previously scanned
/// in the process of plotOn or createHistogram.
/// To do this, you can do the following after the plot has been made:
/// ~~~ {.cpp}
/// profile, RooRealVar * poi, RooCurve * curve ){
///RooCurve * curve = 0;
/// ~~~

TList* ProfileInspector::GetListOfProfilePlots( RooAbsData& data, RooStats::ModelConfig * config)
{

  const RooArgSet* poi_set = config->GetParametersOfInterest();
  const RooArgSet* nuis_params=config->GetNuisanceParameters();
  RooAbsPdf* pdf =  config->GetPdf();


  if(!poi_set){
    cout << "no parameters of interest" << endl;
    return 0;
  }

  if(poi_set->getSize()!=1){
    cout << "only one parameter of interest is supported currently" << endl;
    return 0;
  }
  RooRealVar* poi = (RooRealVar*) poi_set->first();


  if(!nuis_params){
    cout << "no nuisance parameters" << endl;
    return 0;
  }

  if(!pdf){
    cout << "pdf not set" << endl;
    return 0;
  }

  RooAbsReal* nll = pdf->createNLL(data);
  RooAbsReal* profile = nll->createProfile(*poi);

  TList * list = new TList;
  Int_t curve_N=100;
  Double_t* curve_x=0;
//   if(curve){
//     curve_N=curve->GetN();
//     curve_x=curve->GetX();
//     } else {
  Double_t max = dynamic_cast<RooAbsRealLValue*>(poi)->getMax();
  Double_t min = dynamic_cast<RooAbsRealLValue*>(poi)->getMin();
  Double_t step = (max-min)/(curve_N-1);
  curve_x=new Double_t[curve_N];
  for(int i=0; i<curve_N; ++i){
     curve_x[i]=min+step*i;
  }
//   }

  map<string, std::vector<Double_t> > name_val;
  for(int i=0; i<curve_N; i++){
    poi->setVal(curve_x[i]);
    profile->getVal();

    TIterator* nuis_params_itr=nuis_params->createIterator();
    TObject* nuis_params_obj;
    while((nuis_params_obj=nuis_params_itr->Next())){
       RooRealVar* nuis_param = dynamic_cast<RooRealVar*>(nuis_params_obj);
       if(nuis_param) {
          string name = nuis_param->GetName();
          if(nuis_params->getSize()==0) continue;
          if(nuis_param && (! nuis_param->isConstant())){
             if(name_val.find(name)==name_val.end()) name_val[name]=std::vector<Double_t>(curve_N);
             name_val[name][i]=nuis_param->getVal();

             if(i==curve_N-1){
                TGraph* g = new TGraph(curve_N, curve_x, &(name_val[name].front()));
                g->SetName((name+"_"+string(poi->GetName())+"_profile").c_str());
                g->GetXaxis()->SetTitle(poi->GetName());
                g->GetYaxis()->SetTitle(nuis_param->GetName());
                g->SetTitle("");
                list->Add(g);
             }
          }
       }
    }
  }

  delete [] curve_x;


  delete nll;
  delete profile;
  return list;
}
