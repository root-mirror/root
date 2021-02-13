// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodFDA                                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Function discriminant analysis (FDA). This simple classifier              *
 *      fits any user-defined TFormula (via option configuration string) to       *
 *      the training data by requiring a formula response of 1 (0) to signal      *
 *      (background) events. The parameter fitting is done via the abstract       *
 *      class FitterBase, featuring Monte Carlo sampling, Genetic                 *
 *      Algorithm, Simulated Annealing, MINUIT and combinations of these.         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker  <Andreas.Hocker@cern.ch> - CERN, Switzerland             *
 *      Peter Speckmayer <speckmay@mail.cern.ch>  - CERN, Switzerland             *
 *                                                                                *
 * Copyright (c) 2005-2010:                                                       *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodFDA
#define ROOT_TMVA_MethodFDA

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodFDA                                                            //
//                                                                      //
// Function discriminant analysis (FDA). This simple classifier         //
// fits any user-defined TFormula (via option configuration string) to  //
// the training data by requiring a formula response of 1 (0) to signal //
// (background) events. The parameter fitting is done via the abstract  //
// class FitterBase, featuring Monte Carlo sampling, Genetic            //
// Algorithm, Simulated Annealing, MINUIT and combinations of these.    //
//                                                                      //
// Can compute one-dimensional regression                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/MethodBase.h"
#include "TMVA/IFitterTarget.h"
#include <vector>

class TFormula;

namespace TMVA {

   class Interval;
   class Event;
   class FitterBase;

   class MethodFDA : public MethodBase, public IFitterTarget {

   public:

      MethodFDA( const TString& jobName,
                 const TString& methodTitle,
                 DataSetInfo& theData,
                 const TString& theOption = "");

      MethodFDA( DataSetInfo& theData,
                 const TString& theWeightFile);

      virtual ~MethodFDA( void );

      Bool_t HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets );

      // training method
      void Train( void );

      using MethodBase::ReadWeightsFromStream;

      void AddWeightsXMLTo      ( void* parent     ) const;

      void ReadWeightsFromStream( std::istream & i );
      void ReadWeightsFromXML   ( void* wghtnode );

      // calculate the MVA value
      Double_t GetMvaValue( Double_t* err = 0, Double_t* errUpper = 0 );

      virtual const std::vector<Float_t>& GetRegressionValues();
      virtual const std::vector<Float_t>& GetMulticlassValues();

      void Init( void );

      // ranking of input variables
      const Ranking* CreateRanking() { return 0; }

      Double_t EstimatorFunction( std::vector<Double_t>& );

      // no check of options at this place
      void CheckSetup() {}

   protected:

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      void MakeClassSpecific( std::ostream&, const TString& ) const;

      // get help message text
      void GetHelpMessage() const;

   private:

      // compute multiclass values
      void CalculateMulticlassValues( const TMVA::Event*& evt, std::vector<Double_t>& parameters, std::vector<Float_t>& values);


      // create and interpret formula expression and compute estimator
      void     CreateFormula   ();
      Double_t InterpretFormula( const Event*, std::vector<Double_t>::iterator begin, std::vector<Double_t>::iterator end );

      // clean up
      void ClearAll();

      // print fit results
      void PrintResults( const TString&, std::vector<Double_t>&, const Double_t ) const;

      // the option handling methods
      void DeclareOptions();
      void ProcessOptions();

      TString                fFormulaStringP;     // string with function
      TString                fParRangeStringP;    // string with ranges of parameters
      TString                fFormulaStringT;     // string with function
      TString                fParRangeStringT;    // string with ranges of parameters

      TFormula*              fFormula;            // the discrimination function
      UInt_t                 fNPars;              // number of parameters
      std::vector<Interval*> fParRange;           // ranges of parameters
      std::vector<Double_t>  fBestPars;           // the pars that optimise (minimise) the estimator
      TString                fFitMethod;          // estimator optimisation method
      TString                fConverger;          // fitmethod uses fConverger as intermediate step to converge into local minimas
      FitterBase*            fFitter;             // the fitter used in the training
      IFitterTarget*         fConvergerFitter;    // intermediate fitter


      // sum of weights (this should become centrally available through the dataset)
      Double_t               fSumOfWeightsSig;    // sum of weights (signal)
      Double_t               fSumOfWeightsBkg;    // sum of weights (background)
      Double_t               fSumOfWeights;       // sum of weights

      //
      Int_t                  fOutputDimensions;   // number of output values

      ClassDef(MethodFDA,0);  // Function Discriminant Analysis
   };

} // namespace TMVA

#endif // MethodFDA_H
