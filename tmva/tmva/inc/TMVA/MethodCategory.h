// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss,Or Cohen

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodCompositeBase                                                   *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Virtual base class for all MVA method                                     *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker    <Andreas.Hocker@cern.ch> - CERN, Switzerland           *
 *      Joerg Stelzer      <Joerg.Stelzer@cern.ch>  - CERN, Switzerland           *
 *      Peter Speckmayer   <Peter.Speckmayer@cern.ch> - CERN, Switzerland         *
 *      Helge Voss         <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany   *
 *      Eckhard v. Toerne  <evt@uni-bonn.de>        - U of Bonn, Germany          *
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
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodCategory
#define ROOT_TMVA_MethodCategory

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodCategory                                                       //
//                                                                      //
// Class for categorizing the phase space                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iosfwd>
#include <vector>

#include "TMVA/MethodBase.h"

#include "TMVA/MethodCompositeBase.h"

namespace TMVA {

   class Factory;  // DSMTEST
   class Reader;   // DSMTEST
   class MethodBoost;   // DSMTEST
   class DataSetManager;  // DSMTEST
   namespace Experimental {
   class Classification;
   }
   class MethodCategory : public MethodCompositeBase {
      friend class Experimental::Classification;

   public :

      // constructors
      MethodCategory( const TString& jobName,
                      const TString& methodTitle,
                      DataSetInfo& theData,
                      const TString& theOption = "" );

      MethodCategory( DataSetInfo& dsi,
                      const TString& theWeightFile );

      virtual ~MethodCategory( void );

      virtual Bool_t HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ );
      // training and boosting all the classifiers
      void Train( void );

      // ranking of input variables
      const Ranking* CreateRanking();

      // saves the name and options string of the boosted classifier
      TMVA::IMethod* AddMethod(const TCut&,
                               const TString& theVariables,
                               Types::EMVA theMethod,
                               const TString& theTitle,
                               const TString& theOptions);

      void AddWeightsXMLTo( void* parent ) const;
      void ReadWeightsFromXML( void* wghtnode );

      Double_t GetMvaValue( Double_t* err=0, Double_t* errUpper = 0 );

      // regression response
      virtual const std::vector<Float_t>& GetRegressionValues();

      // multi class response
      virtual const std::vector<Float_t> &GetMulticlassValues();

      virtual void MakeClass( const TString& = TString("") ) const {};

   private :

      // initializing mostly monitoring tools of the category process
      void Init();

      // the option handling methods
      void DeclareOptions();
      void ProcessOptions();

      // build the cut formula for event categorization
      Bool_t PassesCut( const Event* ev, UInt_t methodIdx );

   protected:

      // vectors that contain the added methods and the cuts on which they are to be called
      std::vector<IMethod*>               fMethods;
      std::vector<TCut>                   fCategoryCuts;
      std::vector<UInt_t>                 fCategorySpecIdx;
      std::vector<TString>                fVars;
      std::vector <std::vector <UInt_t> > fVarMaps;

      // get help message text
      void GetHelpMessage() const;

      TMVA::DataSetInfo& CreateCategoryDSI(const TCut&, const TString&, const TString&);

   private:

      void InitCircularTree(const DataSetInfo& dsi);

      TTree *                    fCatTree; //! needed in conjunction with TTreeFormulas for evaluation category expressions
      std::vector<TTreeFormula*> fCatFormulas;

      DataSetManager* fDataSetManager; // DSMTEST
      friend class Factory; // DSMTEST
      friend class Reader;  // DSMTEST
      friend class MethodBoost;  // DSMTEST

      ClassDef(MethodCategory,0);
   };
}

#endif
