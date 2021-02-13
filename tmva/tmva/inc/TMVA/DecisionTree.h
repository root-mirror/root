// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss, Jan Therhaag, Eckhard von Toerne 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : DecisionTree                                                          *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation of a Decision Tree                                         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *      Jan Therhaag       <Jan.Therhaag@cern.ch>     - U of Bonn, Germany        *
 *      Eckhard v. Toerne  <evt@uni-bonn.de>          - U of Bonn, Germany        *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         * 
 *      U. of Victoria, Canada                                                    * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://mva.sourceforge.net/license.txt)                                       *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_DecisionTree
#define ROOT_TMVA_DecisionTree

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// DecisionTree                                                         //
//                                                                      //
// Implementation of a Decision Tree                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TH2.h"
#include <vector>

#include "TMVA/Types.h"
#include "TMVA/DecisionTreeNode.h"
#include "TMVA/BinaryTree.h"
#include "TMVA/BinarySearchTree.h"
#include "TMVA/SeparationBase.h"
#include "TMVA/RegressionVariance.h"
#include "TMVA/DataSetInfo.h"

#ifdef R__USE_IMT
#include <ROOT/TThreadExecutor.hxx>
#include "TSystem.h"
#endif

class TRandom3;

namespace TMVA {

   class Event;

   class DecisionTree : public BinaryTree {

   private:

      static const Int_t fgRandomSeed; // set nonzero for debugging and zero for random seeds

   public:

      typedef std::vector<TMVA::Event*> EventList;
      typedef std::vector<const TMVA::Event*> EventConstList;

      // the constructur needed for the "reading" of the decision tree from weight files
      DecisionTree( void );

      // the constructur needed for constructing the decision tree via training with events
      DecisionTree( SeparationBase *sepType, Float_t minSize,
                    Int_t nCuts, DataSetInfo* = NULL,
                    UInt_t cls =0,
                    Bool_t randomisedTree=kFALSE, Int_t useNvars=0, Bool_t usePoissonNvars=kFALSE, 
                    UInt_t nMaxDepth=9999999, 
                    Int_t iSeed=fgRandomSeed, Float_t purityLimit=0.5,
                    Int_t treeID = 0);

      // copy constructor
      DecisionTree (const DecisionTree &d);

      virtual ~DecisionTree( void );

      // Retrieves the address of the root node
      virtual DecisionTreeNode* GetRoot() const { return static_cast<TMVA::DecisionTreeNode*>(fRoot); }
      virtual DecisionTreeNode * CreateNode(UInt_t) const { return new DecisionTreeNode(); }
      virtual BinaryTree* CreateTree() const { return new DecisionTree(); }
      static  DecisionTree* CreateFromXML(void* node, UInt_t tmva_Version_Code = TMVA_VERSION_CODE);
      virtual const char* ClassName() const { return "DecisionTree"; }

      // building of a tree by recursivly splitting the nodes

      //      UInt_t BuildTree( const EventList & eventSample,
      //                        DecisionTreeNode *node = NULL);
      UInt_t BuildTree( const EventConstList & eventSample,
                        DecisionTreeNode *node = NULL);
      // determine the way how a node is split (which variable, which cut value)

      Double_t TrainNode( const EventConstList & eventSample,  DecisionTreeNode *node ) { return TrainNodeFast( eventSample, node ); }
      Double_t TrainNodeFast( const EventConstList & eventSample,  DecisionTreeNode *node );
      Double_t TrainNodeFull( const EventConstList & eventSample,  DecisionTreeNode *node );
      void    GetRandomisedVariables(Bool_t *useVariable, UInt_t *variableMap, UInt_t & nVars);
      std::vector<Double_t>  GetFisherCoefficients(const EventConstList &eventSample, UInt_t nFisherVars, UInt_t *mapVarInFisher);
    
      // fill at tree with a given structure already (just see how many signa/bkgr
      // events end up in each node

      void FillTree( const EventList & eventSample);

      // fill the existing the decision tree structure by filling event
      // in from the top node and see where they happen to end up
      void FillEvent( const TMVA::Event & event,
                      TMVA::DecisionTreeNode *node  );
    
      // returns: 1 = Signal (right),  -1 = Bkg (left)

      Double_t CheckEvent( const TMVA::Event * , Bool_t UseYesNoLeaf = kFALSE ) const;     
      TMVA::DecisionTreeNode* GetEventNode(const TMVA::Event & e) const;

      // return the individual relative variable importance 
      std::vector< Double_t > GetVariableImportance();

      Double_t GetVariableImportance(UInt_t ivar);
    
      // clear the tree nodes (their S/N, Nevents etc), just keep the structure of the tree

      void ClearTree();
    
      // set pruning method
      enum EPruneMethod { kExpectedErrorPruning=0, kCostComplexityPruning, kNoPruning };
      void SetPruneMethod( EPruneMethod m = kCostComplexityPruning ) { fPruneMethod = m; }
    
      // recursive pruning of the tree, validation sample required for automatic pruning
      Double_t PruneTree( const EventConstList* validationSample = NULL );
    
      // manage the pruning strength parameter (iff < 0 -> automate the pruning process)
      void SetPruneStrength( Double_t p ) { fPruneStrength = p; }
      Double_t GetPruneStrength( ) const { return fPruneStrength; }

      // apply pruning validation sample to a decision tree
      void ApplyValidationSample( const EventConstList* validationSample ) const;
    
      // return the misclassification rate of a pruned tree
      Double_t TestPrunedTreeQuality( const DecisionTreeNode* dt = NULL, Int_t mode=0 ) const;
    
      // pass a single validation event throught a pruned decision tree
      void CheckEventWithPrunedTree( const TMVA::Event* ) const;

      // calculate the normalization factor for a pruning validation sample
      Double_t GetSumWeights( const EventConstList* validationSample ) const;
    
      void SetNodePurityLimit( Double_t p ) { fNodePurityLimit = p; }
      Double_t GetNodePurityLimit( ) const { return fNodePurityLimit; }

      void DescendTree( Node *n = NULL );
      void SetParentTreeInNodes( Node *n = NULL );
        
      // retrieve node from the tree. Its position (up to a maximal tree depth of 64)
      // is coded as a sequence of left-right moves starting from the root, coded as
      // 0-1 bit patterns stored in the "long-integer" together with the depth
      Node* GetNode( ULong_t sequence, UInt_t depth );
    
      UInt_t CleanTree(DecisionTreeNode *node=NULL);
     
      void PruneNode(TMVA::DecisionTreeNode *node);    
    
      // prune a node from the tree without deleting its descendants; allows one to
      // effectively prune a tree many times without making deep copies
      void PruneNodeInPlace( TMVA::DecisionTreeNode* node );

      Int_t GetNNodesBeforePruning(){return (fNNodesBeforePruning)?fNNodesBeforePruning:fNNodesBeforePruning=GetNNodes();}
      

      UInt_t CountLeafNodes(TMVA::Node *n = NULL);

      void  SetTreeID(Int_t treeID){fTreeID = treeID;};
      Int_t GetTreeID(){return fTreeID;};

      Bool_t DoRegression() const { return fAnalysisType == Types::kRegression; }
      void SetAnalysisType (Types::EAnalysisType t) { fAnalysisType = t;}
      Types::EAnalysisType GetAnalysisType ( void ) { return fAnalysisType;}
      inline void SetUseFisherCuts(Bool_t t=kTRUE)  { fUseFisherCuts = t;}
      inline void SetMinLinCorrForFisher(Double_t min){fMinLinCorrForFisher = min;}
      inline void SetUseExclusiveVars(Bool_t t=kTRUE){fUseExclusiveVars = t;}
      inline void SetNVars(Int_t n){fNvars = n;}

   private:
      // utility functions
     
      // calculate the Purity out of the number of sig and bkg events collected
      // from individual samples.
    
      // calculates the purity S/(S+B) of a given event sample
      Double_t SamplePurity(EventList eventSample);

      UInt_t    fNvars;          // number of variables used to separate S and B
      Int_t     fNCuts;          // number of grid point in variable cut scans
      Bool_t    fUseFisherCuts;  // use multivariate splits using the Fisher criterium
      Double_t  fMinLinCorrForFisher; // the minimum linear correlation between two variables demanded for use in fisher criterium in node splitting
      Bool_t    fUseExclusiveVars; // individual variables already used in fisher criterium are not anymore analysed individually for node splitting

      SeparationBase *fSepType;  // the separation crition
      RegressionVariance *fRegType;  // the separation crition used in Regression
    
      Double_t  fMinSize;        // min number of events in node
      Double_t  fMinNodeSize;    // min fraction of training events in node
      Double_t  fMinSepGain;     // min number of separation gain to perform node splitting
    
      Bool_t    fUseSearchTree;  // cut scan done with binary trees or simple event loop.
      Double_t  fPruneStrength;  // a parameter to set the "amount" of pruning..needs to be adjusted 
    
      EPruneMethod fPruneMethod; // method used for prunig 
      Int_t    fNNodesBeforePruning; //remember this one (in case of pruning, it allows to monitor the before/after

      Double_t  fNodePurityLimit;// purity limit to decide whether a node is signal
    
      Bool_t    fRandomisedTree; // choose at each node splitting a random set of variables 
      Int_t     fUseNvars;       // the number of variables used in randomised trees;
      Bool_t    fUsePoissonNvars; // use "fUseNvars" not as fixed number but as mean of a possion distr. in each split
    
      TRandom3  *fMyTrandom;     // random number generator for randomised trees
    
      std::vector< Double_t > fVariableImportance; // the relative importance of the different variables 

      UInt_t     fMaxDepth;      // max depth
      UInt_t     fSigClass;      // class which is treated as signal when building the tree
      static const Int_t  fgDebugLevel = 0;     // debug level determining some printout/control plots etc.
      Int_t     fTreeID;        // just an ID number given to the tree.. makes debugging easier as tree knows who he is.

      Types::EAnalysisType  fAnalysisType;   // kClassification(=0=false) or kRegression(=1=true)

      DataSetInfo*  fDataSetInfo;

      ClassDef(DecisionTree,0);               // implementation of a Decision Tree
   };
  
} // namespace TMVA

#endif 
