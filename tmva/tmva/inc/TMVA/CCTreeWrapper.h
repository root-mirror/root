
/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : CCTreeWrapper                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description: a light wrapper of a decision tree, used to perform cost          *
 *              complexity pruning "in-place" Cost Complexity Pruning             *
 *                                                                                *  
 * Author: Doug Schouten (dschoute@sfu.ca)                                        *
 *                                                                                *
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

#ifndef ROOT_TMVA_CCTreeWrapper
#define ROOT_TMVA_CCTreeWrapper

#include "TMVA/Event.h"
#include "TMVA/SeparationBase.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/DataSet.h"
#include "TMVA/Version.h"
#include <vector>
#include <string>
#include <sstream>

namespace TMVA {

   class CCTreeWrapper {

   public:

      typedef std::vector<Event*> EventList;

      /////////////////////////////////////////////////////////////
      // CCTreeNode - a light wrapper of a decision tree node    //
      //                                                         //
      /////////////////////////////////////////////////////////////

      class CCTreeNode : virtual public Node {

      public:

         CCTreeNode( DecisionTreeNode* n = NULL );
         virtual ~CCTreeNode( );
      
         virtual Node* CreateNode() const { return new CCTreeNode(); }

         // set |~T_t|, the number of terminal descendants of node t 
         inline void SetNLeafDaughters( Int_t N ) { fNLeafDaughters = (N > 0 ? N : 0); }

         // return |~T_t|
         inline Int_t GetNLeafDaughters() const { return fNLeafDaughters; }

         // set R(t), the node resubstitution estimate (Gini, misclassification, etc.) for the node t
         inline void SetNodeResubstitutionEstimate( Double_t R ) { fNodeResubstitutionEstimate = (R >= 0 ? R : 0.0); }
      
         // return R(t) for node t
         inline Double_t GetNodeResubstitutionEstimate( ) const { return fNodeResubstitutionEstimate; }

         // set R(T_t) = sum[t' in ~T_t]{ R(t) }, the resubstitution estimate for the branch rooted at
         // node t (it is an estimate because it is calculated from the training dataset, i.e., the original tree)
         inline void SetResubstitutionEstimate( Double_t R ) { fResubstitutionEstimate = (R >= 0 ?  R : 0.0); }
      
         // return R(T_t) for node t
         inline Double_t GetResubstitutionEstimate( ) const { return fResubstitutionEstimate; }
      
         // set the critical point of alpha
         //             R(t) - R(T_t)
         //  alpha_c <  ------------- := g(t)
         //              |~T_t| - 1
         // which is the value of alpha such that the branch rooted at node t is pruned
         inline void SetAlphaC( Double_t alpha ) { fAlphaC = alpha; }

         // get the critical alpha value for this node
         inline Double_t GetAlphaC( ) const { return fAlphaC; }

         // set the minimum critical alpha value for descendants of node t ( G(t) = min(alpha_c, g(t_l), g(t_r)) )
         inline void SetMinAlphaC( Double_t alpha ) { fMinAlphaC = alpha; }

         // get the minimum critical alpha value 
         inline Double_t GetMinAlphaC( ) const { return fMinAlphaC; }

         // get the pointer to the wrapped DT node
         inline DecisionTreeNode* GetDTNode( ) const { return fDTNode; }

         // get pointers to children, mother in the CC tree
         inline CCTreeNode* GetLeftDaughter( ) { return dynamic_cast<CCTreeNode*>(GetLeft()); }
         inline CCTreeNode* GetRightDaughter( ) { return dynamic_cast<CCTreeNode*>(GetRight()); }
         inline CCTreeNode* GetMother( ) { return dynamic_cast<CCTreeNode*>(GetParent()); }

         // printout of the node (can be read in with ReadDataRecord)
         virtual void Print( std::ostream& os ) const;

         // recursive printout of the node and its daughters 
         virtual void PrintRec ( std::ostream& os ) const;

         virtual void AddAttributesToNode(void* node) const;
         virtual void AddContentToNode(std::stringstream& s) const;
         

         // test event if it decends the tree at this node to the right  
         inline virtual Bool_t GoesRight( const Event& e ) const { return (GetDTNode() != NULL ? 
                                                                           GetDTNode()->GoesRight(e) : false); }
      
         // test event if it decends the tree at this node to the left 
         inline virtual Bool_t GoesLeft ( const Event& e ) const { return (GetDTNode() != NULL ? 
                                                                           GetDTNode()->GoesLeft(e) : false); }
         // initialize a node from a data record
         virtual void ReadAttributes(void* node, UInt_t tmva_Version_Code = TMVA_VERSION_CODE);
         virtual void ReadContent(std::stringstream& s);
         virtual Bool_t ReadDataRecord( std::istream& in, UInt_t tmva_Version_Code = TMVA_VERSION_CODE );
      
      private:

         Int_t fNLeafDaughters; //! number of terminal descendants
         Double_t fNodeResubstitutionEstimate; //! R(t) = misclassification rate for node t
         Double_t fResubstitutionEstimate; //! R(T_t) = sum[t' in ~T_t]{ R(t) }
         Double_t fAlphaC; //! critical point, g(t) = alpha_c(t)
         Double_t fMinAlphaC; //! G(t), minimum critical point of t and its descendants
         DecisionTreeNode* fDTNode; //! pointer to wrapped node in the decision tree
      };

      CCTreeWrapper( DecisionTree* T,  SeparationBase* qualityIndex );
      ~CCTreeWrapper( );

      // return the decision tree output for an event 
      Double_t CheckEvent( const TMVA::Event & e, Bool_t useYesNoLeaf = false );
      // return the misclassification rate of a pruned tree for a validation event sample
      Double_t TestTreeQuality( const EventList* validationSample );
      Double_t TestTreeQuality( const DataSet* validationSample );

      // remove the branch rooted at node t
      void PruneNode( CCTreeNode* t );
      // initialize the node t and all its descendants
      void InitTree( CCTreeNode* t );

      // return the root node for this tree
      CCTreeNode* GetRoot() { return fRoot; }
   private:
      SeparationBase* fQualityIndex;  //! pointer to the used quality index calculator
      DecisionTree* fDTParent;        //! pointer to underlying DecisionTree
      CCTreeNode* fRoot;              //! the root node of the (wrapped) decision Tree
   };

}

#endif



