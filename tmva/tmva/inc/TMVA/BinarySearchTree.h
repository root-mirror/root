// @(#)root/tmva $Id$    
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : BinarySearchTree                                                      *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      BinarySearchTree incl. volume Search method                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
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
 **********************************************************************************/

#ifndef ROOT_TMVA_BinarySearchTree
#define ROOT_TMVA_BinarySearchTree

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// BinarySearchTree                                                     //
//                                                                      //
// A simple Binary search tree including volume search method           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <queue>
#include <utility>

#ifndef ROOT_time
#include "time.h"
#endif

#include "TMVA/Volume.h"
#include "TMVA/BinaryTree.h"
#include "TMVA/BinarySearchTreeNode.h"
#include "TMVA/VariableInfo.h"

class TString;
class TTree;

// -----------------------------------------------------------------------------
// the binary search tree

namespace TMVA {

   class Event;
   //   class MethodBase;
   
   class BinarySearchTree : public BinaryTree {
      
   public:
      
      // constructor
      BinarySearchTree( void );
    
      // copy constructor
      BinarySearchTree (const BinarySearchTree &b);

      // destructor
      virtual ~BinarySearchTree( void );
    
      virtual Node * CreateNode( UInt_t ) const { return new BinarySearchTreeNode(); }
      virtual BinaryTree* CreateTree() const { return new BinarySearchTree(); }
      static BinarySearchTree* CreateFromXML(void* node, UInt_t tmva_Version_Code = TMVA_VERSION_CODE);
      virtual const char* ClassName() const { return "BinarySearchTree"; }

      // Searches for a node with the specified data 
      // by calling  the private, recursive, function for searching
      BinarySearchTreeNode* Search( Event * event ) const;
    
      // Adds an item to the tree, 
      void Insert( const Event * );
    
      // get sum of weights of the nodes;
      Double_t GetSumOfWeights( void ) const;

      //get sum of weights of the nodes of given type;
      Double_t GetSumOfWeights( Int_t theType ) const;
    
      //set the periode (number of variables)
      void SetPeriode( Int_t p )      { fPeriod = p; }

      // return periode (number of variables)
      UInt_t  GetPeriode( void ) const { return fPeriod; }

      // counts events (weights) within a given volume 
      Double_t SearchVolume( Volume*, std::vector<const TMVA::BinarySearchTreeNode*>* events = 0 );
    
      // Create the search tree from the event collection 
      // using ONLY the variables specified in "theVars"
      Double_t Fill( const std::vector<TMVA::Event*>& events, const std::vector<Int_t>& theVars, Int_t theType = -1 );
    
      // create the search tree from the events in a TTree
      // using ALL the variables specified included in the Event
      Double_t Fill( const std::vector<TMVA::Event*>& events, Int_t theType = -1 );

      void NormalizeTree ();      
      
      void CalcStatistics( TMVA::Node* n = 0 );      
      void Clear         ( TMVA::Node* n = 0 );

      // access to mean for signal and background for each variable
      Float_t Mean(Types::ESBType sb, UInt_t var ) { return fMeans[sb==Types::kSignal?0:1][var]; }

      // access to RMS for signal and background for each variable
      Float_t RMS(Types::ESBType sb, UInt_t var ) { return fRMS[sb==Types::kSignal?0:1][var]; }

      // access to Minimum for signal and background for each variable
      Float_t Min(Types::ESBType sb, UInt_t var ) { return fMin[sb==Types::kSignal?0:1][var]; }

      // access to Maximum for signal and background for each variable
      Float_t Max(Types::ESBType sb, UInt_t var ) { return fMax[sb==Types::kSignal?0:1][var]; }

      Int_t SearchVolumeWithMaxLimit( TMVA::Volume*, std::vector<const TMVA::BinarySearchTreeNode*>* events = 0, Int_t = -1);

      // access to RMS for each variable
      Float_t RMS(UInt_t var ) { return fRMS[0][var]; } // attention! class 0 is taken as signal!

      void SetNormalize( Bool_t norm ) { fCanNormalize = norm; }

   private:

      // add a new  node to the tree (as daughter) 
      void       Insert( const Event*, Node* );
      // recursively search the nodes for Event
      BinarySearchTreeNode*      Search( Event*, Node *) const ;
    
      //check of Event variables lie with the volumde
      Bool_t   InVolume    (const std::vector<Float_t>&, Volume* ) const;
      //
      void     DestroyNode ( BinarySearchTreeNode* );


      void     NormalizeTree( std::vector< std::pair< Double_t, const TMVA::Event* > >::iterator, 
                              std::vector< std::pair< Double_t, const TMVA::Event* > >::iterator, UInt_t );

      // recursive search through daughter nodes in weight counting
      Double_t SearchVolume( Node*, Volume*, Int_t, 
                             std::vector<const TMVA::BinarySearchTreeNode*>* events );
      UInt_t fPeriod;            // periode (number of event variables)
      UInt_t fCurrentDepth;      // internal variable, counting the depth of the tree during insertion    
      Bool_t fStatisticsIsValid; // flag if last stat calculation is still valid, set to false if new node is insert

      std::vector<Float_t>        fMeans[2];    // mean for signal and background for each variable
      std::vector<Float_t>        fRMS[2];      // RMS for signal and background for each variable
      std::vector<Float_t>        fMin[2];      // RMS for signal and background for each variable
      std::vector<Float_t>        fMax[2];      // RMS for signal and background for each variable
      std::vector<Double_t>       fSum[2];      // Sum for signal and background for each variable
      std::vector<Double_t>       fSumSq[2];    // Squared Sum for signal and background for each variable
      Double_t                    fNEventsW[2]; // Number of events per class, taking into account event weights
      Double_t                    fSumOfWeights;// Total number of events (weigthed) counted during filling
                                                // should be the same as fNEventsW[0]+fNEventsW[1].. used as a check
      
      Bool_t                      fCanNormalize; // the tree can be normalised
      std::vector< std::pair<Double_t,const TMVA::Event*> > fNormalizeTreeTable;
      
      ClassDef(BinarySearchTree,0); // Binary search tree including volume search method  
   };
  
} // namespace TMVA

#endif
