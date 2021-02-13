// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Classes: Node                                                                  *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Node for the BinarySearch or Decision Trees                               *
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
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_Node
#define ROOT_TMVA_Node

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Node                                                                 //
//                                                                      //
// Node base class for the BinarySearch or Decision Trees Nodes         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iosfwd>
#include <string>
#include <sstream>

#include "Rtypes.h"
#include "TMVA/Version.h"

namespace TMVA {

   class Node;
   class Event;
   class BinaryTree;

   std::ostream& operator<<( std::ostream& os, const Node& node );
   std::ostream& operator<<( std::ostream& os, const Node* node );

   // a class used to identify a Node; (needed for recursive reading from text file)
   // (currently it is NOT UNIQUE... but could eventually made it
   // a node in the tree structure
   class Node {

      // output operator for a node
      friend std::ostream& operator << (std::ostream& os, const Node& node);
      // output operator with a pointer to the node (which still prints the node itself)
      friend std::ostream& operator << (std::ostream& os, const Node* node);

   public:

      // constructor of a node
      Node();

      // constructor of a daughter node as a daughter of 'p'
      Node( Node* p, char pos );

      // copy constructor
      Node( const Node &n );

      // destructor
      virtual ~Node();

      virtual Node* CreateNode() const = 0;

      // test event if i{ descends the tree at this node to the right
      virtual Bool_t GoesRight( const Event& ) const = 0;
      // test event if it descends the tree at this node to the left

      virtual Bool_t GoesLeft ( const Event& ) const = 0;
      // test event if it is equal to the event that "makes the node" (just for the "search tree"

      // return pointer to the left/right daughter or parent node
      inline virtual Node* GetLeft  () const { return fLeft;   }
      inline virtual Node* GetRight () const { return fRight;  }
      inline virtual Node* GetParent() const { return fParent; }

      // set pointer to the left/right daughter or parent node
      inline virtual void SetLeft  (Node* l) { fLeft   = l;}
      inline virtual void SetRight (Node* r) { fRight  = r;}
      inline virtual void SetParent(Node* p) { fParent = p;}

      //recursively go through the part of the tree below this node and count all daughters
      Int_t  CountMeAndAllDaughters() const;

      // printout of the node
      virtual void Print( std::ostream& os ) const = 0;

      // recursive printout of the node and it daughters
      virtual void PrintRec ( std::ostream& os ) const = 0;

      void* AddXMLTo(void* parent) const;
      void  ReadXML(void* node, UInt_t tmva_Version_Code = TMVA_VERSION_CODE );
      virtual void AddAttributesToNode(void* node) const = 0;
      virtual void AddContentToNode(std::stringstream& s) const = 0;

      // Set depth, layer of the where the node is within the tree, seen from the top (root)
      void SetDepth(UInt_t d){fDepth=d;}

      // Return depth, layer of the where the node is within the tree, seen from the top (root)
      UInt_t GetDepth() const {return fDepth;}

      // set node position, i.e, the node is a left (l) or right (r) daughter
      void SetPos(char s) {fPos=s;}

      // Return the node position, i.e, the node is a left (l) or right (r) daughter
      char GetPos() const {return fPos;}

      // Return the pointer to the Parent tree to which the Node belongs
      virtual TMVA::BinaryTree* GetParentTree() const {return fParentTree;}

      // set the pointer to the Parent Tree to which the Node belongs
      virtual void SetParentTree(TMVA::BinaryTree* t) {fParentTree = t;}

      int GetCount();

      virtual Bool_t ReadDataRecord( std::istream&, UInt_t tmva_Version_Code = TMVA_VERSION_CODE ) = 0;
      virtual void ReadAttributes(void* node, UInt_t tmva_Version_Code = TMVA_VERSION_CODE  ) = 0;
      virtual void ReadContent(std::stringstream& s) =0;

   protected:

      Node*   fParent;              // the previous (parent) node
      Node*   fLeft;                // pointers to the two "daughter" nodes
      Node*   fRight;               // pointers to the two "daughter" nodes

      char    fPos;                 // position, i.e. it is a left (l) or right (r) daughter
      UInt_t  fDepth;               // depth of the node within the tree (seen from root node)

      BinaryTree*  fParentTree;     // pointer to the parent tree to which the Node belongs
   private:

      static Int_t fgCount;         // counter of all nodes present.. for debug.. to spot memory leaks...

   public:
      ClassDef(Node,0); // Node for the BinarySearch or Decision Trees
   };

} // namespace TMVA

#endif

