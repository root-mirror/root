/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::DecisionTree                                                    *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      IPruneTool - a helper interface class to prune a decision tree            *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Doug Schouten <dschoute@sfu.ca> - Simon Fraser U., Canada                 *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://mva.sourceforge.net/license.txt)                                       *
 **********************************************************************************/

#ifndef ROOT_TMVA_IPruneTool
#define ROOT_TMVA_IPruneTool

#include <iosfwd>
#include <vector>

#include "TMVA/SeparationBase.h"

#include "TMVA/DecisionTree.h"

namespace TMVA {

   //    class MsgLogger;

   ////////////////////////////////////////////////////////////
   // Basic struct for saving relevant pruning information   //
   ////////////////////////////////////////////////////////////
   class PruningInfo {

   public:

   PruningInfo( ) : QualityIndex(0), PruneStrength(0), PruneSequence(0) {}
      PruningInfo( Double_t q, Double_t alpha, std::vector<DecisionTreeNode*> sequence );
      Double_t QualityIndex; //! quality measure for a pruned subtree T of T_max
      Double_t PruneStrength; //! the regularization parameter for pruning
      std::vector<DecisionTreeNode*> PruneSequence; //! the sequence of pruning locations in T_max that yields T
   };

   inline PruningInfo::PruningInfo( Double_t q, Double_t alpha, std::vector<DecisionTreeNode*> sequence )
      : QualityIndex(q), PruneStrength(alpha), PruneSequence(sequence) {}

/*! \class TMVA::IPruneTool
\ingroup TMVA
IPruneTool - a helper interface class to prune a decision tree

Any tool which implements the interface should provide two modes for tree pruning:

  1. automatically find the "best" prune strength by minimizing the error rate on a test sample
     if SetAutomatic() is called, or if automatic = kTRUE argument is set in CalculatePruningInfo()
     In this case, the PruningInfo object returned contains the error rate of the optimally pruned
     tree, the optimal prune strength, and the sequence of nodes to prune to obtain the optimal
     pruned tree from the original DecisionTree

  2. a user-provided pruning strength parameter is used to prune the tree, in which case the returned
     PruningInfo object has QualityIndex = -1, PruneStrength = user prune strength, and PruneSequence
     is the list of nodes to prune
*/

   class IPruneTool {

   public:

      typedef std::vector<const Event*> EventSample;

      IPruneTool( );
      virtual ~IPruneTool();

   public:

      // returns the PruningInfo object for a given tree and test sample
      virtual PruningInfo* CalculatePruningInfo( DecisionTree* dt, const EventSample* testEvents = NULL,
                                                 Bool_t isAutomatic = kFALSE ) = 0;

   public:

      // set the prune strength parameter to use in pruning
      inline void SetPruneStrength( Double_t alpha ) { fPruneStrength = alpha; }
      // return the prune strength the tool is using
      inline Double_t GetPruneStrength( ) const { return fPruneStrength; }

      // if the prune strength parameter is < 0, the tool will automatically find an optimal strength
      // set the tool to automatic mode
      inline void SetAutomatic( ) { fPruneStrength = -1.0; };
      inline Bool_t IsAutomatic( ) const { return fPruneStrength <= 0.0; }

   protected:

      //       mutable MsgLogger* fLogger; //! output stream to save logging information
      //       MsgLogger& Log() const { return *fLogger; }
      Double_t fPruneStrength; //! regularization parameter in pruning


      Double_t S, B;
   };

   inline IPruneTool::IPruneTool( ) :
      fPruneStrength(0.0),
      S(0),
      B(0)
         {}
}

#endif
