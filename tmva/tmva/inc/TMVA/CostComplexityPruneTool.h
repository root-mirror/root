/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::DecisionTree                                                    *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation of a Decision Tree                                         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *      Doug Schouten   <dschoute@sfu.ca>        - Simon Fraser U., Canada        *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://mva.sourceforge.net/license.txt)                                       *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_CostComplexityPruneTool
#define ROOT_TMVA_CostComplexityPruneTool

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CostComplexityPruneTool - a class to prune a decision tree using the Cost Complexity method            //
// (see "Classification and Regression Trees" by Leo Breiman et al)                                       //
//                                                                                                        //
// Some definitions:                                                                                      //
//                                                                                                        //
// T_max - the initial, usually highly overtrained tree, that is to be pruned back                        //
// R(T) - quality index (Gini, misclassification rate, or other) of a tree T                              //
// ~T - set of terminal nodes in T                                                                        //
// T' - the pruned subtree of T_max that has the best quality index R(T')                                 //
// alpha - the prune strength parameter in Cost Complexity pruning (R_alpha(T) = R(T) + alpha*|~T|)       //
//                                                                                                        //
// There are two running modes in CostComplexityPruneTool: (i) one may select a prune strength and prune  //
// the tree T_max until the criterion                                                                     //
//             R(T) - R(t)                                                                                //
//  alpha <    ----------                                                                                 //
//             |~T_t| - 1                                                                                 //
//                                                                                                        //
// is true for all nodes t in T, or (ii) the algorithm finds the sequence of critical points              //
// alpha_k < alpha_k+1 ... < alpha_K such that T_K = root(T_max) and then selects the optimally-pruned    //
// subtree, defined to be the subtree with the best quality index for the validation sample.              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TMVA/SeparationBase.h"
#include "TMVA/GiniIndex.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/Event.h"
#include "TMVA/IPruneTool.h"
#include <vector>

namespace TMVA {

   class CostComplexityPruneTool : public IPruneTool {
   public:
      CostComplexityPruneTool( SeparationBase* qualityIndex = NULL );
      virtual ~CostComplexityPruneTool( );

      // calculate the prune sequence for a given tree
      virtual PruningInfo* CalculatePruningInfo( DecisionTree* dt, const IPruneTool::EventSample* testEvents = NULL, Bool_t isAutomatic = kFALSE );

   private:
      SeparationBase* fQualityIndexTool; //! the quality index used to calculate R(t), R(T) = sum[t in ~T]{ R(t) }

      std::vector<DecisionTreeNode*> fPruneSequence; //! map of weakest links (i.e., branches to prune) -> pruning index
      std::vector<Double_t> fPruneStrengthList; //! map of alpha -> pruning index
      std::vector<Double_t> fQualityIndexList; //! map of R(T) -> pruning index

      Int_t fOptimalK; //! the optimal index of the prune sequence

   private:
      // set the meta data used for cost complexity pruning
      void InitTreePruningMetaData( DecisionTreeNode* n );

      // optimize the pruning sequence
      void Optimize( DecisionTree* dt, Double_t weights );

      mutable MsgLogger* fLogger; //! output stream to save logging information
      MsgLogger& Log() const { return *fLogger; }

   };
}


#endif
