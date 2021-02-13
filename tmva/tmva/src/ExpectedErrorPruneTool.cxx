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

/*! \class TMVA::ExpectedErrorPruneTool
\ingroup TMVA

A helper class to prune a decision tree using the expected error (C4.5) method

Uses an upper limit on the error made by the classification done by each node.
If the \f$ \frac{S}{S+B} \f$ of the node is \f$ f \f$, then according to the
training sample, the error rate (fraction of misclassified events by this
node) is \f$ (1-f) \f$. Now \f$ f \f$ has a statistical error according to the
binomial distribution hence the error on \f$ f \f$ can be estimated (same error
as the binomial error for efficiency calculations
\f$ (\sigma = \sqrt{\frac{(eff(1-eff)}{nEvts}}) \f$

This tool prunes branches from a tree if the expected error of a node is less
than that of the sum of the error in its descendants.

*/

#include "TMVA/ExpectedErrorPruneTool.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/IPruneTool.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Types.h"

#include "RtypesCore.h"
#include "Rtypes.h"
#include "TMath.h"

#include <map>

// pin the vtable here.
TMVA::IPruneTool::~IPruneTool() {}

////////////////////////////////////////////////////////////////////////////////

TMVA::ExpectedErrorPruneTool::ExpectedErrorPruneTool() :
   IPruneTool(),
   fDeltaPruneStrength(0),
   fNodePurityLimit(1),
   fLogger( new MsgLogger("ExpectedErrorPruneTool") )
{}

////////////////////////////////////////////////////////////////////////////////

TMVA::ExpectedErrorPruneTool::~ExpectedErrorPruneTool()
{
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////

TMVA::PruningInfo*
TMVA::ExpectedErrorPruneTool::CalculatePruningInfo( DecisionTree* dt,
                                                    const IPruneTool::EventSample* validationSample,
                                                    Bool_t isAutomatic )
{
   if( isAutomatic ) {
      //SetAutomatic( );
      isAutomatic = kFALSE;
      Log() << kWARNING << "Sorry automatic pruning strength determination is not implemented yet" << Endl;
   }
   if( dt == NULL || (IsAutomatic() && validationSample == NULL) ) {
      // must have a valid decision tree to prune, and if the prune strength
      // is to be chosen automatically, must have a test sample from
      // which to calculate the quality of the pruned tree(s)
      return NULL;
   }
   fNodePurityLimit = dt->GetNodePurityLimit();

   if(IsAutomatic()) {
      Log() << kFATAL << "Sorry automatic pruning strength determination is not implemented yet" << Endl;
      /*
        dt->ApplyValidationSample(validationSample);
        Double_t weights = dt->GetSumWeights(validationSample);
        // set the initial prune strength
        fPruneStrength = 1.0e-3; //! FIXME somehow make this automatic, it depends strongly on the tree structure
        // better to set it too small, it will be increased automatically
        fDeltaPruneStrength = 1.0e-5;
        Int_t nnodes = this->CountNodes((DecisionTreeNode*)dt->GetRoot());

        Bool_t forceStop = kFALSE;
        Int_t errCount = 0,
        lastNodeCount = nnodes;

        // find the maximum prune strength that still leaves the root's daughter nodes

        while ( nnodes > 1 && !forceStop ) {
        fPruneStrength += fDeltaPruneStrength;
        Log() << "----------------------------------------------------" << Endl;
        FindListOfNodes((DecisionTreeNode*)dt->GetRoot());
        for( UInt_t i = 0; i < fPruneSequence.size(); i++ )
        fPruneSequence[i]->SetTerminal(); // prune all the nodes from the sequence
        // test the quality of the pruned tree
        Double_t quality = 1.0 - dt->TestPrunedTreeQuality()/weights;
        fQualityMap.insert(std::make_pair<const Double_t,Double_t>(quality,fPruneStrength));

        nnodes = CountNodes((DecisionTreeNode*)dt->GetRoot()); // count the number of nodes in the pruned tree

        Log() << "Prune strength : " << fPruneStrength << Endl;
        Log() << "Had " << lastNodeCount << " nodes, now have " << nnodes << Endl;
        Log() << "Quality index is: " << quality << Endl;

        if (lastNodeCount == nnodes) errCount++;
        else {
        errCount=0; // reset counter
        if ( nnodes < lastNodeCount / 2 ) {
        Log() << "Decreasing fDeltaPruneStrength to " << fDeltaPruneStrength/2.0
        << " because the number of nodes in the tree decreased by a factor of 2." << Endl;
        fDeltaPruneStrength /= 2.;
        }
        }
        lastNodeCount = nnodes;
        if (errCount > 20) {
        Log() << "Increasing fDeltaPruneStrength to " << fDeltaPruneStrength*2.0
        << " because the number of nodes in the tree didn't change." << Endl;
        fDeltaPruneStrength *= 2.0;
        }
        if (errCount > 40) {
        Log() << "Having difficulty determining the optimal prune strength, bailing out!" << Endl;
        forceStop = kTRUE;
        }
        // reset the tree for the next iteration
        for( UInt_t i = 0; i < fPruneSequence.size(); i++ )
        fPruneSequence[i]->SetTerminal(false);
        fPruneSequence.clear();
        }
        // from the set of pruned trees, find the one with the optimal quality index
        std::multimap<Double_t,Double_t>::reverse_iterator it = fQualityMap.rend(); ++it;
        fPruneStrength = it->second;
        FindListOfNodes((DecisionTreeNode*)dt->GetRoot());

        // adjust the step size for the next tree automatically
        fPruneStrength = 1.0e-3;
        fDeltaPruneStrength = (fPruneStrength - 1.0)/(Double_t)fQualityMap.size();

        return new PruningInfo(it->first, it->second, fPruneSequence);
      */
      return NULL;
   }
   else { // no automatic pruning - just use the provided prune strength parameter
      FindListOfNodes( (DecisionTreeNode*)dt->GetRoot() );
      return new PruningInfo( -1.0, fPruneStrength, fPruneSequence );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// recursive pruning of nodes using the Expected Error Pruning (EEP)

void TMVA::ExpectedErrorPruneTool::FindListOfNodes( DecisionTreeNode* node )
{
   TMVA::DecisionTreeNode *l = (TMVA::DecisionTreeNode*)node->GetLeft();
   TMVA::DecisionTreeNode *r = (TMVA::DecisionTreeNode*)node->GetRight();
   if (node->GetNodeType() == 0 && !(node->IsTerminal())) { // check all internal nodes
      this->FindListOfNodes(l);
      this->FindListOfNodes(r);
      if (this->GetSubTreeError(node) >= this->GetNodeError(node)) {
         //node->Print(Log());
         fPruneSequence.push_back(node);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// calculate the expected statistical error on the subtree below "node"
/// which is used in the expected error pruning

Double_t TMVA::ExpectedErrorPruneTool::GetSubTreeError( DecisionTreeNode* node ) const
{
   DecisionTreeNode *l = (DecisionTreeNode*)node->GetLeft();
   DecisionTreeNode *r = (DecisionTreeNode*)node->GetRight();
   if (node->GetNodeType() == 0 && !(node->IsTerminal())) {
      Double_t subTreeError =
         (l->GetNEvents() * this->GetSubTreeError(l) +
          r->GetNEvents() * this->GetSubTreeError(r)) /
         node->GetNEvents();
      return subTreeError;
   }
   else {
      return this->GetNodeError(node);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate an UPPER limit on the error made by the classification done
/// by this node. If the S/S+B of the node is f, then according to the
/// training sample, the error rate (fraction of misclassified events by
/// this node) is (1-f)
/// Now f has a statistical error according to the binomial distribution
/// hence the error on f can be estimated (same error as the binomial error
/// for efficiency calculations
/// \f$ (\sigma = \sqrt{\frac{(eff(1-eff)}{nEvts}}) \f$

Double_t TMVA::ExpectedErrorPruneTool::GetNodeError( DecisionTreeNode *node ) const
{
   Double_t errorRate = 0;

   Double_t nEvts = node->GetNEvents();

   // fraction of correctly classified events by this node:
   Double_t f = 0;
   if (node->GetPurity() > fNodePurityLimit) f = node->GetPurity();
   else  f = (1-node->GetPurity());

   Double_t df = TMath::Sqrt(f*(1-f)/nEvts);

   errorRate = std::min(1.0,(1.0 - (f-fPruneStrength*df)));

   // -------------------------------------------------------------------
   // standard algorithm:
   // step 1: Estimate error on node using Laplace estimate
   //         NodeError = (N - n + k -1 ) / (N + k)
   //   N: number of events
   //   k: number of event classes (2 for Signal, Background)
   //   n: n event out of N belong to the class which has the majority in the node
   // step 2: Approximate "backed-up" error assuming we did not prune
   //   (I'm never quite sure if they consider whole subtrees, or only 'next-to-leaf'
   //    nodes)...
   //   Subtree error = Sum_children ( P_i * NodeError_i)
   //    P_i = probability of the node to make the decision, i.e. fraction of events in
   //          leaf node ( N_leaf / N_parent)
   // step 3:

   // Minimum Error Pruning (MEP) according to Niblett/Bratko
   //# of correctly classified events by this node:
   //Double_t n=f*nEvts ;
   //Double_t p_apriori = 0.5, m=100;
   //errorRate = (nEvts  - n + (1-p_apriori) * m ) / (nEvts  + m);

   // Pessimistic error Pruning (proposed by Quinlan (error estimat with continuity approximation)
   //# of correctly classified events by this node:
   //Double_t n=f*nEvts ;
   //errorRate = (nEvts  - n + 0.5) / nEvts ;

   //const Double Z=.65;
   //# of correctly classified events by this node:
   //Double_t n=f*nEvts ;
   //errorRate = (f + Z*Z/(2*nEvts ) + Z*sqrt(f/nEvts  - f*f/nEvts  + Z*Z/4/nEvts /nEvts ) ) / (1 + Z*Z/nEvts );
   //errorRate = (n + Z*Z/2 + Z*sqrt(n - n*n/nEvts  + Z*Z/4) )/ (nEvts  + Z*Z);
   //errorRate = 1 - errorRate;
   // -------------------------------------------------------------------

   return errorRate;
}

