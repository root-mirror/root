#ifndef ROOT_TMVA_CCPruner
#define ROOT_TMVA_CCPruner
/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : CCPruner                                                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description: Cost Complexity Pruning                                           *
 * 
 * Author: Doug Schouten (dschoute@sfu.ca)
 *
 *                                                                                *
 * Copyright (c) 2007:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Texas at Austin, USA                                                *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CCPruner - a helper class to prune a decision tree using the Cost Complexity method                    //
// (see Classification and Regression Trees by Leo Breiman et al)                                         //
//                                                                                                        //
// Some definitions:                                                                                      //
//                                                                                                        //
// T_max - the initial, usually highly overtrained tree, that is to be pruned back                        // 
// R(T) - quality index (Gini, misclassification rate, or other) of a tree T                              //
// ~T - set of terminal nodes in T                                                                        //
// T' - the pruned subtree of T_max that has the best quality index R(T')                                 //
// alpha - the prune strength parameter in Cost Complexity pruning (R_alpha(T) = R(T) + alpha// |~T|)     //
//                                                                                                        //
// There are two running modes in CCPruner: (i) one may select a prune strength and prune back            //
// the tree T_max until the criterion                                                                     //
//             R(T) - R(t)                                                                                //
//  alpha <    ----------                                                                                 //
//             |~T_t| - 1                                                                                 //
//                                                                                                        //
// is true for all nodes t in T, or (ii) the algorithm finds the sequence of critical points              //
// alpha_k < alpha_k+1 ... < alpha_K such that T_K = root(T_max) and then selects the optimally-pruned    //
// subtree, defined to be the subtree with the best quality index for the validation sample.              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "TMVA/DecisionTree.h"

/* #ifndef ROOT_TMVA_DecisionTreeNode */
/* #include "TMVA/DecisionTreeNode.h" */
/* #endif */

#include "TMVA/Event.h"
#include <vector>

namespace TMVA {
   class DataSet;
   class DecisionTreeNode;
   class SeparationBase;

   class CCPruner {
   public: 
      typedef std::vector<Event*> EventList;

      CCPruner( DecisionTree* t_max, 
                const EventList* validationSample,
                SeparationBase* qualityIndex = NULL );

      CCPruner( DecisionTree* t_max, 
                const DataSet* validationSample,
                SeparationBase* qualityIndex = NULL );

      ~CCPruner( );

      // set the pruning strength parameter alpha (if alpha < 0, the optimal alpha is calculated)
      void SetPruneStrength( Float_t alpha = -1.0 );

      void Optimize( );

      // return the list of pruning locations to define the optimal subtree T' of T_max
      std::vector<TMVA::DecisionTreeNode*> GetOptimalPruneSequence( ) const; 

      // return the quality index from the validation sample for the optimal subtree T'
      inline Float_t GetOptimalQualityIndex( ) const { return (fOptimalK >= 0 && fQualityIndexList.size() > 0 ?
                                                               fQualityIndexList[fOptimalK] : -1.0); }

      // return the prune strength (=alpha) corresponding to the prune sequence
      inline Float_t GetOptimalPruneStrength( ) const { return (fOptimalK >= 0 && fPruneStrengthList.size() > 0 ?
                                                                fPruneStrengthList[fOptimalK] : -1.0); }
   
   private:
      Float_t              fAlpha; //! regularization parameter in CC pruning
      const EventList*     fValidationSample; //! the event sample to select the optimally-pruned tree
      const DataSet*       fValidationDataSet; //! the event sample to select the optimally-pruned tree
      SeparationBase*      fQualityIndex; //! the quality index used to calculate R(t), R(T) = sum[t in ~T]{ R(t) }
      Bool_t               fOwnQIndex; //! flag indicates if fQualityIndex is owned by this

      DecisionTree*        fTree; //! (pruned) decision tree

      std::vector<TMVA::DecisionTreeNode*> fPruneSequence; //! map of weakest links (i.e., branches to prune) -> pruning index
      std::vector<Float_t> fPruneStrengthList;  //! map of alpha -> pruning index
      std::vector<Float_t> fQualityIndexList;   //! map of R(T) -> pruning index

      Int_t                fOptimalK;           //! index of the optimal tree in the pruned tree sequence
      Bool_t               fDebug;              //! debug flag
   };
}

inline void TMVA::CCPruner::SetPruneStrength( Float_t alpha ) {
   fAlpha = (alpha > 0 ? alpha : 0.0);
}
    

#endif


