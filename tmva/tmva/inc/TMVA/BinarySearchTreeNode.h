// @(#)root/tmva $Id$    
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Classes: Node, NodeID                                                          *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Node for the BinarySearch                                                 *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
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

#ifndef ROOT_TMVA_BinarySearchTreeNode
#define ROOT_TMVA_BinarySearchTreeNode

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// BinarySearchTreeNode                                                 //
//                                                                      //
// Node for the BinarySearch  Tree                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iosfwd>
#include <vector>
#include <string>
#include <sstream>
#include "Rtypes.h"

#include "TMVA/Node.h"

namespace TMVA {

   class Event;

   // a class used to identify a Node; (needed for recursive reading from text file)
   // (currently it is NOT UNIQUE... but could eventually made it
   // a node in the tree structure
   class BinarySearchTreeNode : public Node  {

   public:

      // constructor of a node for the search tree
      BinarySearchTreeNode( const Event* e = NULL, UInt_t signalClass=0 );

      // constructor of a daughter node as a daughter of 'p'
      BinarySearchTreeNode( BinarySearchTreeNode* parent, char pos );

      // copy constructor
      BinarySearchTreeNode ( const BinarySearchTreeNode &n, 
                             BinarySearchTreeNode* parent = NULL);

      // destructor
      virtual ~BinarySearchTreeNode ();

      virtual Node* CreateNode() const { return new BinarySearchTreeNode(); }

      // test event if it decends the tree at this node to the right  
      virtual Bool_t GoesRight( const Event& ) const;
      // test event if it decends the tree at this node to the left 

      virtual Bool_t GoesLeft ( const Event& ) const;
      // test event if it is equal to the event that "makes the node" (just for the "search tree"  

      virtual Bool_t EqualsMe ( const Event& ) const;

      // set index of variable used for discrimination at this node
      inline void SetSelector( Short_t i) { fSelector = i; }
      // return index of variable used for discrimination at this node 
      inline Short_t GetSelector() const { return fSelector; }

      const std::vector<Float_t> & GetEventV() const { return fEventV; }
      Float_t                      GetWeight() const { return fWeight; }
      UInt_t                       GetClass()  const { return fClass; }
      //      Bool_t                       IsSignal()  const { return (fClass == fSignalClass); }

      const std::vector<Float_t> & GetTargets() const { return fTargets; }


      // printout of the node
      virtual void Print( std::ostream& os ) const;

      // recursive printout of the node and it daughters 
      virtual void PrintRec( std::ostream& os ) const;

      virtual void AddAttributesToNode(void* node) const;
      virtual void AddContentToNode(std::stringstream& s) const;

      // Read the data block
      virtual Bool_t ReadDataRecord( std::istream& is, UInt_t tmva_Version_Code = TMVA_VERSION_CODE );
      virtual void ReadAttributes(void* node, UInt_t tmva_Version_Code = TMVA_VERSION_CODE );
      virtual void ReadContent(std::stringstream& s);

   private: 
      std::vector<Float_t> fEventV;
      std::vector<Float_t> fTargets;

      Float_t     fWeight;
      UInt_t      fClass;

      Short_t     fSelector;       // index of variable used in node selection (decision tree) 

      ClassDef(BinarySearchTreeNode,0); // Node for the BinarySearchTree
   };

} // namespace TMVA

#endif
