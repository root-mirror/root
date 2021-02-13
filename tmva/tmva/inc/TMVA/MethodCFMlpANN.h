// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodCFMlpANN                                                        *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Interface for Clermond-Ferrand artificial neural network.                 *
 *      The ANN code has been translated from FORTRAN77 (f2c);                    *
 *      see files: MethodCFMlpANN_f2c_mlpl3.cpp                                   *
 *                 MethodCFMlpANN_f2c_datacc.cpp                                  *
 *                                                                                *
 *      --------------------------------------------------------------------      *
 *      Reference for the original FORTRAN version:                               *
 *           Authors  : J. Proriol and contributions from ALEPH-Clermont-Fd       *
 *                      Team members. Contact : gaypas@afal11.cern.ch             *
 *                                                                                *
 *           Copyright: Laboratoire Physique Corpusculaire                        *
 *                      Universite de Blaise Pascal, IN2P3/CNRS                   *
 *      --------------------------------------------------------------------      *
 *                                                                                *
 * Usage: options are given through Factory:                                      *
 *            factory->BookMethod( "MethodCFMlpANN", OptionsString );             *
 *                                                                                *
 *        where:                                                                  *
 *            TString OptionsString = "n_training_cycles:n_hidden_layers"         *
 *                                                                                *
 *        default is:  n_training_cycles = 5000, n_layers = 4                     *
 *        note that the number of hidden layers in the NN is                      *
 *                                                                                *
 *            n_hidden_layers = n_layers - 2                                      *
 *                                                                                *
 *        since there is one input and one output layer. The number of            *
 *        nodes (neurons) is predefined to be                                     *
 *                                                                                *
 *           n_nodes[i] = nvars + 1 - i (where i=1..n_layers)                     *
 *                                                                                *
 *        with nvars being the number of variables used in the NN.                *
 *        Hence, the default case is: n_neurons(layer 1 (input)) : nvars          *
 *                                    n_neurons(layer 2 (hidden)): nvars-1        *
 *                                    n_neurons(layer 3 (hidden)): nvars-1        *
 *                                    n_neurons(layer 4 (out))   : 2              *
 *                                                                                *
 *        This artificial neural network usually needs a relatively large         *
 *        number of cycles to converge (8000 and more). Overtraining can          *
 *        be efficiently tested by comparing the signal and background            *
 *        output of the NN for the events that were used for training and         *
 *        an independent data sample (with equal properties). If the separation   *
 *        performance is significantly better for the training sample, the        *
 *        NN interprets statistical effects, and is hence overtrained. In         *
 *        this case, the number of cycles should be reduced, or the size          *
 *        of the training sample increased.                                       *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Xavier Prudent  <prudent@lapp.in2p3.fr>  - LAPP, France                   *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodCFMlpANN
#define ROOT_TMVA_MethodCFMlpANN

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodCFMlpANN                                                       //
//                                                                      //
// Interface for Clermond-Ferrand artificial neural network             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iosfwd>
#include <vector>

#include "TMVA/MethodBase.h"
#include "TMVA/MethodCFMlpANN_Utils.h"
#include "TMatrixF.h"

namespace TMVA {

   class MethodCFMlpANN : public MethodBase, MethodCFMlpANN_Utils {

   public:

      MethodCFMlpANN( const TString& jobName,
                      const TString& methodTitle,
                      DataSetInfo& theData,
                      const TString& theOption = "3000:N-1:N-2");

      MethodCFMlpANN( DataSetInfo& theData,
                      const TString& theWeightFile);

      virtual ~MethodCFMlpANN( void );

      virtual Bool_t HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ );

      // training method
      void Train( void );

      using MethodBase::ReadWeightsFromStream;

      // write weights to file
      void AddWeightsXMLTo( void* parent ) const;

      // read weights from file
      void ReadWeightsFromStream( std::istream& istr );
      void ReadWeightsFromXML( void* wghtnode );
      // calculate the MVA value
      Double_t GetMvaValue( Double_t* err = 0, Double_t* errUpper = 0 );

      // data accessors for external functions
      Double_t GetData ( Int_t isel, Int_t ivar ) const { return (*fData)(isel, ivar); }
      Int_t    GetClass( Int_t ivar             ) const { return (*fClass)[ivar]; }


      // ranking of input variables
      const Ranking* CreateRanking() { return 0; }

   protected:

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      void MakeClassSpecific( std::ostream&, const TString& ) const;

      // header and auxiliary classes
      void MakeClassSpecificHeader( std::ostream&, const TString& = "" ) const;

      // get help message text
      void GetHelpMessage() const;

      Int_t DataInterface( Double_t*, Double_t*, Int_t*, Int_t*, Int_t*, Int_t*,
                           Double_t*, Int_t*, Int_t* );

   private:

      void PrintWeights( std::ostream & o ) const;

      // the option handling methods
      void DeclareOptions();
      void ProcessOptions();

      // LUTs
      TMatrixF       *fData;     // the (data,var) string
      std::vector<Int_t> *fClass;    // the event class (1=signal, 2=background)

      Int_t         fNlayers;   // number of layers (including input and output layers)
      Int_t         fNcycles;   // number of training cycles
      Int_t*        fNodes;     // number of nodes per layer

      // additional member variables for the independent NN::Evaluation phase
      Double_t**    fYNN;       // weights
      TString       fLayerSpec; // the hidden layer specification string
      Int_t MethodCFMlpANN_nsel;

      // auxiliary member functions
      Double_t EvalANN( std::vector<Double_t>&, Bool_t& isOK );
      void     NN_ava ( Double_t* );
      Double_t NN_fonc( Int_t, Double_t ) const;

      // default initialisation
      void Init( void );

      ClassDef(MethodCFMlpANN,0); // Interface for Clermond-Ferrand artificial neural network
   };

} // namespace TMVA

#endif
