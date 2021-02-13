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

#ifndef ROOT_TMVA_ExpectedErrorPruneTool
#define ROOT_TMVA_ExpectedErrorPruneTool

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ExpectedErrorPruneTool - a helper class to prune a decision tree using the expected error (C4.5) method //
//                                                                                                         //
// Uses an upper limit on the error made by the classification done by each node. If the S/S+B of the node //
// is f, then according to the training sample, the error rate (fraction of misclassified events by this   //
// node) is (1-f). Now f has a statistical error according to the binomial distribution hence the error on //
// f can be estimated (same error as the binomial error for efficiency calculations                         //
// ( sigma = sqrt(eff(1-eff)/nEvts ) )                                                                     //
//                                                                                                         //
// This tool prunes branches from a tree if the expected error of a node is less than that of the sum  of  //
// the error in its descendants.                                                                           //
//                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "TMath.h"

#include "TMVA/IPruneTool.h"

namespace TMVA {

   class MsgLogger;

   class ExpectedErrorPruneTool : public IPruneTool {
   public:
      ExpectedErrorPruneTool( );
      virtual ~ExpectedErrorPruneTool( );

      // returns the PruningInfo object for a given tree and test sample
      virtual PruningInfo* CalculatePruningInfo( DecisionTree* dt, const IPruneTool::EventSample* testEvents = NULL,
                                                 Bool_t isAutomatic = kFALSE );

   public:
      // set the increment dalpha with which to scan for the optimal prune strength
      inline void SetPruneStrengthIncrement( Double_t dalpha ) { fDeltaPruneStrength = dalpha; }

   private:
      void FindListOfNodes( DecisionTreeNode* node );
      Double_t GetNodeError( DecisionTreeNode* node ) const;
      Double_t GetSubTreeError( DecisionTreeNode* node ) const;
      Int_t CountNodes( DecisionTreeNode* node, Int_t icount = 0 );

      Double_t fDeltaPruneStrength; //! the stepsize for optimizing the pruning strength parameter
      Double_t fNodePurityLimit; //! the purity limit for labelling a terminal node as signal
      std::vector<DecisionTreeNode*> fPruneSequence; //! the (optimal) prune sequence
      //      std::multimap<const Double_t, Double_t> fQualityMap; //! map of tree quality <=> prune strength
      mutable MsgLogger* fLogger;   // message logger
      MsgLogger& Log() const { return *fLogger; }
   };

   inline Int_t ExpectedErrorPruneTool::CountNodes( DecisionTreeNode* node, Int_t icount ) {
      DecisionTreeNode* l = (DecisionTreeNode*)node->GetLeft();
      DecisionTreeNode* r = (DecisionTreeNode*)node->GetRight();
      Int_t counter = icount + 1; // count this node
      if(!(node->IsTerminal()) && l != NULL && r != NULL) {
         counter = CountNodes(l,counter);
         counter = CountNodes(r,counter);
      }
      return counter;
   }
}

#endif

