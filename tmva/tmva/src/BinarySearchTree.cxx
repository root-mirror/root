// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : BinarySearchTree                                                      *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header file for description)                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <stelzer@cern.ch>        - DESY, Germany                  *
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
 *                                                                                *
 **********************************************************************************/

/*! \class TMVA::BinarySearchTree
\ingroup TMVA

A simple Binary search tree including a volume search method.

*/

#include <stdexcept>
#include <cstdlib>
#include <queue>
#include <algorithm>

#include "TMath.h"

#include "TMatrixDBase.h"
#include "TTree.h"

#include "TMVA/MsgLogger.h"
#include "TMVA/MethodBase.h"
#include "TMVA/Tools.h"
#include "TMVA/Event.h"
#include "TMVA/BinarySearchTree.h"

#include "TMVA/BinaryTree.h"
#include "TMVA/Types.h"
#include "TMVA/Node.h"

ClassImp(TMVA::BinarySearchTree);

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TMVA::BinarySearchTree::BinarySearchTree( void ) :
BinaryTree(),
   fPeriod      ( 1 ),
   fCurrentDepth( 0 ),
   fStatisticsIsValid( kFALSE ),
   fSumOfWeights( 0 ),
   fCanNormalize( kFALSE )
{
   fNEventsW[0]=fNEventsW[1]=0.;
}

////////////////////////////////////////////////////////////////////////////////
/// copy constructor that creates a true copy, i.e. a completely independent tree

TMVA::BinarySearchTree::BinarySearchTree( const BinarySearchTree &b)
   : BinaryTree(),
     fPeriod      ( b.fPeriod ),
     fCurrentDepth( 0 ),
     fStatisticsIsValid( kFALSE ),
     fSumOfWeights( b.fSumOfWeights ),
     fCanNormalize( kFALSE )
{
   fNEventsW[0]=fNEventsW[1]=0.;
   Log() << kFATAL << " Copy constructor not implemented yet " << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::BinarySearchTree::~BinarySearchTree( void )
{
   for(std::vector< std::pair<Double_t, const TMVA::Event*> >::iterator pIt = fNormalizeTreeTable.begin();
       pIt != fNormalizeTreeTable.end(); ++pIt) {
      delete pIt->second;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// re-create a new tree (decision tree or search tree) from XML

TMVA::BinarySearchTree* TMVA::BinarySearchTree::CreateFromXML(void* node, UInt_t tmva_Version_Code ) {
   std::string type("");
   gTools().ReadAttr(node,"type", type);
   BinarySearchTree* bt = new BinarySearchTree();
   bt->ReadXML( node, tmva_Version_Code );
   return bt;
}

////////////////////////////////////////////////////////////////////////////////
/// insert a new "event" in the binary tree

void TMVA::BinarySearchTree::Insert( const Event* event )
{
   fCurrentDepth=0;
   fStatisticsIsValid = kFALSE;

   if (this->GetRoot() == NULL) {           // If the list is empty...
      this->SetRoot( new BinarySearchTreeNode(event)); //Make the new node the root.
      // have to use "s" for start as "r" for "root" would be the same as "r" for "right"
      this->GetRoot()->SetPos('s');
      this->GetRoot()->SetDepth(0);
      fNNodes = 1;
      fSumOfWeights = event->GetWeight();
      ((BinarySearchTreeNode*)this->GetRoot())->SetSelector((UInt_t)0);
      this->SetPeriode(event->GetNVariables());
   }
   else {
      // sanity check:
      if (event->GetNVariables() != (UInt_t)this->GetPeriode()) {
         Log() << kFATAL << "<Insert> event vector length != Periode specified in Binary Tree" << Endl
               << "--- event size: " << event->GetNVariables() << " Periode: " << this->GetPeriode() << Endl
               << "--- and all this when trying filling the "<<fNNodes+1<<"th Node" << Endl;
      }
      // insert a new node at the propper position
      this->Insert(event, this->GetRoot());
   }

   // normalise the tree to speed up searches
   if (fCanNormalize) fNormalizeTreeTable.push_back( std::make_pair(0.0,new const Event(*event)) );
}

////////////////////////////////////////////////////////////////////////////////
/// private internal function to insert a event (node) at the proper position

void TMVA::BinarySearchTree::Insert( const Event *event,
                                     Node *node )
{
   fCurrentDepth++;
   fStatisticsIsValid = kFALSE;

   if (node->GoesLeft(*event)){    // If the adding item is less than the current node's data...
      if (node->GetLeft() != NULL){            // If there is a left node...
         // Add the new event to the left node
         this->Insert(event, node->GetLeft());
      }
      else {                             // If there is not a left node...
         // Make the new node for the new event
         BinarySearchTreeNode* current = new BinarySearchTreeNode(event);
         fNNodes++;
         fSumOfWeights += event->GetWeight();
         current->SetSelector(fCurrentDepth%((Int_t)event->GetNVariables()));
         current->SetParent(node);          // Set the new node's previous node.
         current->SetPos('l');
         current->SetDepth( node->GetDepth() + 1 );
         node->SetLeft(current);            // Make it the left node of the current one.
      }
   }
   else if (node->GoesRight(*event)) { // If the adding item is less than or equal to the current node's data...
      if (node->GetRight() != NULL) {              // If there is a right node...
         // Add the new node to it.
         this->Insert(event, node->GetRight());
      }
      else {                                 // If there is not a right node...
         // Make the new node.
         BinarySearchTreeNode* current = new BinarySearchTreeNode(event);
         fNNodes++;
         fSumOfWeights += event->GetWeight();
         current->SetSelector(fCurrentDepth%((Int_t)event->GetNVariables()));
         current->SetParent(node);              // Set the new node's previous node.
         current->SetPos('r');
         current->SetDepth( node->GetDepth() + 1 );
         node->SetRight(current);               // Make it the left node of the current one.
      }
   }
   else Log() << kFATAL << "<Insert> neither left nor right :)" << Endl;
}

////////////////////////////////////////////////////////////////////////////////
///search the tree to find the node matching "event"

TMVA::BinarySearchTreeNode* TMVA::BinarySearchTree::Search( Event* event ) const
{
   return this->Search( event, this->GetRoot() );
}

////////////////////////////////////////////////////////////////////////////////
/// Private, recursive, function for searching.

TMVA::BinarySearchTreeNode* TMVA::BinarySearchTree::Search(Event* event, Node* node) const
{
   if (node != NULL) {               // If the node is not NULL...
      // If we have found the node...
      if (((BinarySearchTreeNode*)(node))->EqualsMe(*event))
         return (BinarySearchTreeNode*)node;                  // Return it
      if (node->GoesLeft(*event))      // If the node's data is greater than the search item...
         return this->Search(event, node->GetLeft());  //Search the left node.
      else                          //If the node's data is less than the search item...
         return this->Search(event, node->GetRight()); //Search the right node.
   }
   else return NULL; //If the node is NULL, return NULL.
}

////////////////////////////////////////////////////////////////////////////////
/// return the sum of event (node) weights

Double_t TMVA::BinarySearchTree::GetSumOfWeights( void ) const
{
   if (fSumOfWeights <= 0) {
      Log() << kWARNING << "you asked for the SumOfWeights, which is not filled yet"
            << " I call CalcStatistics which hopefully fixes things"
            << Endl;
   }
   if (fSumOfWeights <= 0) Log() << kFATAL << " Zero events in your Search Tree" <<Endl;

   return fSumOfWeights;
}

////////////////////////////////////////////////////////////////////////////////
/// return the sum of event (node) weights

Double_t TMVA::BinarySearchTree::GetSumOfWeights( Int_t theType ) const
{
   if (fSumOfWeights <= 0) {
      Log() << kWARNING << "you asked for the SumOfWeights, which is not filled yet"
            << " I call CalcStatistics which hopefully fixes things"
            << Endl;
   }
   if (fSumOfWeights <= 0) Log() << kFATAL << " Zero events in your Search Tree" <<Endl;

   return fNEventsW[ ( theType == Types::kSignal) ? 0 : 1  ];
}

////////////////////////////////////////////////////////////////////////////////
/// create the search tree from the event collection
/// using ONLY the variables specified in "theVars"

Double_t TMVA::BinarySearchTree::Fill( const std::vector<Event*>& events, const std::vector<Int_t>& theVars,
                                       Int_t theType )
{
   fPeriod = theVars.size();
   return Fill(events, theType);
}

////////////////////////////////////////////////////////////////////////////////
/// create the search tree from the events in a TTree
/// using ALL the variables specified included in the Event

Double_t TMVA::BinarySearchTree::Fill( const std::vector<Event*>& events, Int_t theType )
{
   UInt_t n=events.size();

   UInt_t nevents = 0;
   if (fSumOfWeights != 0) {
      Log() << kWARNING
            << "You are filling a search three that is not empty.. "
            << " do you know what you are doing?"
            << Endl;
   }
   for (UInt_t ievt=0; ievt<n; ievt++) {
      // insert event into binary tree
      if (theType == -1 || (Int_t(events[ievt]->GetClass()) == theType) ) {
         this->Insert( events[ievt] );
         nevents++;
         fSumOfWeights += events[ievt]->GetWeight();
      }
   } // end of event loop
   CalcStatistics(0);

   return fSumOfWeights;
}

////////////////////////////////////////////////////////////////////////////////
/// normalises the binary-search tree to reduce the branch length and hence speed up the
/// search procedure (on average).

void TMVA::BinarySearchTree::NormalizeTree ( std::vector< std::pair<Double_t, const TMVA::Event*> >::iterator leftBound,
                                             std::vector< std::pair<Double_t, const TMVA::Event*> >::iterator rightBound,
                                             UInt_t actDim )
{

   if (leftBound == rightBound) return;

   if (actDim == fPeriod)  actDim = 0;
   for (std::vector< std::pair<Double_t, const TMVA::Event*> >::iterator i=leftBound; i!=rightBound; ++i) {
      i->first = i->second->GetValue( actDim );
   }

   std::sort( leftBound, rightBound );

   std::vector< std::pair<Double_t, const TMVA::Event*> >::iterator leftTemp  = leftBound;
   std::vector< std::pair<Double_t, const TMVA::Event*> >::iterator rightTemp = rightBound;

   // meet in the middle
   while (true) {
      --rightTemp;
      if (rightTemp == leftTemp ) {
         break;
      }
      ++leftTemp;
      if (leftTemp  == rightTemp) {
         break;
      }
   }

   std::vector< std::pair<Double_t, const TMVA::Event*> >::iterator mid     = leftTemp;
   std::vector< std::pair<Double_t, const TMVA::Event*> >::iterator midTemp = mid;

   if (mid!=leftBound)--midTemp;

   while (mid != leftBound && mid->second->GetValue( actDim ) == midTemp->second->GetValue( actDim ))  {
      --mid;
      --midTemp;
   }

   Insert( mid->second );

   //    Print(std::cout);
   //    std::cout << std::endl << std::endl;

   NormalizeTree( leftBound, mid, actDim+1 );
   ++mid;
   //    Print(std::cout);
   //    std::cout << std::endl << std::endl;
   NormalizeTree( mid, rightBound, actDim+1 );


   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Normalisation of tree

void TMVA::BinarySearchTree::NormalizeTree()
{
   SetNormalize( kFALSE );
   Clear( NULL );
   this->SetRoot(NULL);
   NormalizeTree( fNormalizeTreeTable.begin(), fNormalizeTreeTable.end(), 0 );
}

////////////////////////////////////////////////////////////////////////////////
/// clear nodes

void TMVA::BinarySearchTree::Clear( Node* n )
{
   BinarySearchTreeNode* currentNode = (BinarySearchTreeNode*)(n == NULL ? this->GetRoot() : n);

   if (currentNode->GetLeft()  != 0) Clear( currentNode->GetLeft()  );
   if (currentNode->GetRight() != 0) Clear( currentNode->GetRight() );

   if (n != NULL) delete n;

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// search the whole tree and add up all weights of events that
/// lie within the given volume

Double_t TMVA::BinarySearchTree::SearchVolume( Volume* volume,
                                               std::vector<const BinarySearchTreeNode*>* events )
{
   return SearchVolume( this->GetRoot(), volume, 0, events );
}

////////////////////////////////////////////////////////////////////////////////
/// recursively walk through the daughter nodes and add up all weights of events that
/// lie within the given volume

Double_t TMVA::BinarySearchTree::SearchVolume( Node* t, Volume* volume, Int_t depth,
                                               std::vector<const BinarySearchTreeNode*>* events )
{
   if (t==NULL) return 0;  // Are we at an outer leave?

   BinarySearchTreeNode* st = (BinarySearchTreeNode*)t;

   Double_t count = 0.0;
   if (InVolume( st->GetEventV(), volume )) {
      count += st->GetWeight();
      if (NULL != events) events->push_back( st );
   }
   if (st->GetLeft()==NULL && st->GetRight()==NULL) {

      return count;  // Are we at an outer leave?
   }

   Bool_t tl, tr;
   Int_t  d = depth%this->GetPeriode();
   if (d != st->GetSelector()) {
      Log() << kFATAL << "<SearchVolume> selector in Searchvolume "
            << d << " != " << "node "<< st->GetSelector() << Endl;
   }
   tl = (*(volume->fLower))[d] <  st->GetEventV()[d];  // Should we descend left?
   tr = (*(volume->fUpper))[d] >= st->GetEventV()[d];  // Should we descend right?

   if (tl) count += SearchVolume( st->GetLeft(),  volume, (depth+1), events );
   if (tr) count += SearchVolume( st->GetRight(), volume, (depth+1), events );

   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// test if the data points are in the given volume

Bool_t TMVA::BinarySearchTree::InVolume(const std::vector<Float_t>& event, Volume* volume ) const
{

   Bool_t result = false;
   for (UInt_t ivar=0; ivar< fPeriod; ivar++) {
      result = ( (*(volume->fLower))[ivar] <  event[ivar] &&
                 (*(volume->fUpper))[ivar] >= event[ivar] );
      if (!result) break;
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// calculate basic statistics (mean, rms for each variable)

void TMVA::BinarySearchTree::CalcStatistics( Node* n )
{
   if (fStatisticsIsValid) return;

   BinarySearchTreeNode * currentNode = (BinarySearchTreeNode*)n;

   // default, start at the tree top, then descend recursively
   if (n == NULL) {
      fSumOfWeights = 0;
      for (Int_t sb=0; sb<2; sb++) {
         fNEventsW[sb]  = 0;
         fMeans[sb]     = std::vector<Float_t>(fPeriod);
         fRMS[sb]       = std::vector<Float_t>(fPeriod);
         fMin[sb]       = std::vector<Float_t>(fPeriod);
         fMax[sb]       = std::vector<Float_t>(fPeriod);
         fSum[sb]       = std::vector<Double_t>(fPeriod);
         fSumSq[sb]     = std::vector<Double_t>(fPeriod);
         for (UInt_t j=0; j<fPeriod; j++) {
            fMeans[sb][j] = fRMS[sb][j] = fSum[sb][j] = fSumSq[sb][j] = 0;
            fMin[sb][j] =  FLT_MAX;
            fMax[sb][j] = -FLT_MAX;
         }
      }
      currentNode = (BinarySearchTreeNode*) this->GetRoot();
      if (currentNode == NULL) return; // no root-node
   }

   const std::vector<Float_t> & evtVec = currentNode->GetEventV();
   Double_t                     weight = currentNode->GetWeight();
   //    Int_t                        type   = currentNode->IsSignal();
   //   Int_t                        type   = currentNode->IsSignal() ? 0 : 1;
   Int_t                        type   = Int_t(currentNode->GetClass())== Types::kSignal ? 0 : 1;

   fNEventsW[type] += weight;
   fSumOfWeights   += weight;

   for (UInt_t j=0; j<fPeriod; j++) {
      Float_t val = evtVec[j];
      fSum[type][j]   += val*weight;
      fSumSq[type][j] += val*val*weight;
      if (val < fMin[type][j]) fMin[type][j] = val;
      if (val > fMax[type][j]) fMax[type][j] = val;
   }

   if ( (currentNode->GetLeft()  != NULL) ) CalcStatistics( currentNode->GetLeft() );
   if ( (currentNode->GetRight() != NULL) ) CalcStatistics( currentNode->GetRight() );

   if (n == NULL) { // i.e. the root node
      for (Int_t sb=0; sb<2; sb++) {
         for (UInt_t j=0; j<fPeriod; j++) {
            if (fNEventsW[sb] == 0) { fMeans[sb][j] = fRMS[sb][j] = 0; continue; }
            fMeans[sb][j] = fSum[sb][j]/fNEventsW[sb];
            fRMS[sb][j]   = TMath::Sqrt(fSumSq[sb][j]/fNEventsW[sb] - fMeans[sb][j]*fMeans[sb][j]);
         }
      }
      fStatisticsIsValid = kTRUE;
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// recursively walk through the daughter nodes and add up all weights of events that
/// lie within the given volume a maximum number of events can be given

Int_t TMVA::BinarySearchTree::SearchVolumeWithMaxLimit( Volume *volume, std::vector<const BinarySearchTreeNode*>* events,
                                                        Int_t max_points )
{
   if (this->GetRoot() == NULL) return 0;  // Are we at an outer leave?

   std::queue< std::pair< const BinarySearchTreeNode*, Int_t > > queue;
   std::pair< const BinarySearchTreeNode*, Int_t > st = std::make_pair( (const BinarySearchTreeNode*)this->GetRoot(), 0 );
   queue.push( st );

   Int_t count = 0;

   while ( !queue.empty() ) {
      st = queue.front(); queue.pop();

      if (count == max_points)
         return count;

      if (InVolume( st.first->GetEventV(), volume )) {
         count++;
         if (NULL != events) events->push_back( st.first );
      }

      Bool_t tl, tr;
      Int_t d = st.second;
      if ( d == Int_t(this->GetPeriode()) ) d = 0;

      if (d != st.first->GetSelector()) {
         Log() << kFATAL << "<SearchVolume> selector in Searchvolume "
               << d << " != " << "node "<< st.first->GetSelector() << Endl;
      }

      tl = (*(volume->fLower))[d] <  st.first->GetEventV()[d] && st.first->GetLeft()  != NULL;  // Should we descend left?
      tr = (*(volume->fUpper))[d] >= st.first->GetEventV()[d] && st.first->GetRight() != NULL;  // Should we descend right?

      if (tl) queue.push( std::make_pair( (const BinarySearchTreeNode*)st.first->GetLeft(), d+1 ) );
      if (tr) queue.push( std::make_pair( (const BinarySearchTreeNode*)st.first->GetRight(), d+1 ) );
   }

   return count;
}
