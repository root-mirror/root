// @(#)root/tmva $Id$
// Author: Rustem Ospanov

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Node                                                                  *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      kd-tree (binary tree) template                                            *
 *                                                                                *
 * Author:                                                                        *
 *      Rustem Ospanov <rustem@fnal.gov> - U. of Texas at Austin, USA             *
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

#ifndef ROOT_TMVA_NodekNN
#define ROOT_TMVA_NodekNN

// C++
#include <cassert>
#include <list>
#include <string>
#include <iostream>
#include <utility>

// ROOT
#include "RtypesCore.h"

/*! \class TMVA::kNN::Node
\ingroup TMVA
This file contains binary tree and global function template
that searches tree for k-nearest neigbors

Node class template parameter T has to provide these functions:
  rtype GetVar(UInt_t) const;
  - rtype is any type convertible to Float_t
  UInt_t GetNVar(void) const;
  rtype GetWeight(void) const;
  - rtype is any type convertible to Double_t

Find function template parameter T has to provide these functions:
(in addition to above requirements)
  rtype GetDist(Float_t, UInt_t) const;
  - rtype is any type convertible to Float_t
  rtype GetDist(const T &) const;
  - rtype is any type convertible to Float_t

  where T::GetDist(Float_t, UInt_t) <= T::GetDist(const T &)
  for any pair of events and any variable number for these events
*/

namespace TMVA
{
   namespace kNN
   {
      template <class T>
         class Node
         {

         public:

            Node(const Node *parent, const T &event, Int_t mod);
            ~Node();

            const Node* Add(const T &event, UInt_t depth);

            void SetNodeL(Node *node);
            void SetNodeR(Node *node);

            const T& GetEvent() const;

            const Node* GetNodeL() const;
            const Node* GetNodeR() const;
            const Node* GetNodeP() const;

            Double_t GetWeight() const;

            Float_t GetVarDis() const;
            Float_t GetVarMin() const;
            Float_t GetVarMax() const;

            UInt_t GetMod() const;

            void Print() const;
            void Print(std::ostream& os, const std::string &offset = "") const;

         private:

            // these methods are private and not implemented by design
            // use provided public constructor for all uses of this template class
            Node();
            Node(const Node &);
            const Node& operator=(const Node &);

         private:

            const Node* fNodeP;

            Node* fNodeL;
            Node* fNodeR;

            const T fEvent;

            const Float_t fVarDis;

            Float_t fVarMin;
            Float_t fVarMax;

            const UInt_t fMod;
         };

      // recursive search for k-nearest neighbor: k = nfind
      template<class T>
         UInt_t Find(std::list<std::pair<const Node<T> *, Float_t> > &nlist,
                     const Node<T> *node, const T &event, UInt_t nfind);

      // recursive search for k-nearest neighbor
      // find k events with sum of event weights >= nfind
      template<class T>
         UInt_t Find(std::list<std::pair<const Node<T> *, Float_t> > &nlist,
                     const Node<T> *node, const T &event, Double_t nfind, Double_t ncurr);

      // recursively travel upward until root node is reached
      template <class T>
         UInt_t Depth(const Node<T> *node);

      // prInt_t node content and content of its children
      //template <class T>
      //std::ostream& operator<<(std::ostream& os, const Node<T> &node);

      //
      // Inlined functions for Node template
      //
      template <class T>
         inline void Node<T>::SetNodeL(Node<T> *node)
         {
            fNodeL = node;
         }

      template <class T>
         inline void Node<T>::SetNodeR(Node<T> *node)
         {
            fNodeR = node;
         }

      template <class T>
         inline const T& Node<T>::GetEvent() const
         {
            return fEvent;
         }

      template <class T>
         inline const Node<T>* Node<T>::GetNodeL() const
         {
            return fNodeL;
         }

      template <class T>
         inline const Node<T>* Node<T>::GetNodeR() const
         {
            return fNodeR;
         }

      template <class T>
         inline const Node<T>* Node<T>::GetNodeP() const
         {
            return fNodeP;
         }

      template <class T>
         inline Double_t Node<T>::GetWeight() const
         {
            return fEvent.GetWeight();
         }

      template <class T>
         inline Float_t Node<T>::GetVarDis() const
         {
            return fVarDis;
         }

      template <class T>
         inline Float_t Node<T>::GetVarMin() const
         {
            return fVarMin;
         }

      template <class T>
         inline Float_t Node<T>::GetVarMax() const
         {
            return fVarMax;
         }

      template <class T>
         inline UInt_t Node<T>::GetMod() const
         {
            return fMod;
         }

      //
      // Inlined global function(s)
      //
      template <class T>
         inline UInt_t Depth(const Node<T> *node)
         {
            if (!node) return 0;
            else return Depth(node->GetNodeP()) + 1;
         }

   } // end of kNN namespace
} // end of TMVA namespace

////////////////////////////////////////////////////////////////////////////////
template<class T>
TMVA::kNN::Node<T>::Node(const Node<T> *parent, const T &event, const Int_t mod)
:fNodeP(parent),
   fNodeL(0),
   fNodeR(0),
   fEvent(event),
   fVarDis(event.GetVar(mod)),
   fVarMin(fVarDis),
   fVarMax(fVarDis),
   fMod(mod)
{}

////////////////////////////////////////////////////////////////////////////////
template<class T>
TMVA::kNN::Node<T>::~Node()
{
   if (fNodeL) delete fNodeL;
   if (fNodeR) delete fNodeR;
}

////////////////////////////////////////////////////////////////////////////////
/// This is Node member function that adds a new node to a binary tree.
/// each node contains maximum and minimum values of splitting variable
/// left or right nodes are added based on value of splitting variable

template<class T>
const TMVA::kNN::Node<T>* TMVA::kNN::Node<T>::Add(const T &event, const UInt_t depth)
{

   assert(fMod == depth % event.GetNVar() && "Wrong recursive depth in Node<>::Add");

   const Float_t value = event.GetVar(fMod);

   fVarMin = std::min(fVarMin, value);
   fVarMax = std::max(fVarMax, value);

   Node<T> *node = 0;
   if (value < fVarDis) {
      if (fNodeL)
         {
            return fNodeL->Add(event, depth + 1);
         }
      else {
         fNodeL = new Node<T>(this, event, (depth + 1) % event.GetNVar());
         node = fNodeL;
      }
   }
   else {
      if (fNodeR) {
         return fNodeR->Add(event, depth + 1);
      }
      else {
         fNodeR = new Node<T>(this, event, (depth + 1) % event.GetNVar());
         node = fNodeR;
      }
   }

   return node;
}

////////////////////////////////////////////////////////////////////////////////
template<class T>
void TMVA::kNN::Node<T>::Print() const
{
   Print(std::cout);
}

////////////////////////////////////////////////////////////////////////////////
template<class T>
void TMVA::kNN::Node<T>::Print(std::ostream& os, const std::string &offset) const
{
   os << offset << "-----------------------------------------------------------" << std::endl;
   os << offset << "Node: mod " << fMod
      << " at " << fVarDis
      << " with weight: " << GetWeight() << std::endl
      << offset << fEvent;

   if (fNodeL) {
      os << offset << "Has left node " << std::endl;
   }
   if (fNodeR) {
      os << offset << "Has right node" << std::endl;
   }

   if (fNodeL) {
      os << offset << "PrInt_t left node " << std::endl;
      fNodeL->Print(os, offset + " ");
   }
   if (fNodeR) {
      os << offset << "PrInt_t right node" << std::endl;
      fNodeR->Print(os, offset + " ");
   }

   if (!fNodeL && !fNodeR) {
      os << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This is a global templated function that searches for k-nearest neighbors.
/// list contains k or less nodes that are closest to event.
/// only nodes with positive weights are added to list.
/// each node contains maximum and minimum values of splitting variable
/// for all its children - this range is checked to avoid descending into
/// nodes that are definitely outside current minimum neighbourhood.
///
/// This function should be modified with care.

template<class T>
UInt_t TMVA::kNN::Find(std::list<std::pair<const TMVA::kNN::Node<T> *, Float_t> > &nlist,
                       const TMVA::kNN::Node<T> *node, const T &event, const UInt_t nfind)
{
   if (!node || nfind < 1) {
      return 0;
   }

   const Float_t value = event.GetVar(node->GetMod());

   if (node->GetWeight() > 0.0) {

      Float_t max_dist = 0.0;

      if (!nlist.empty()) {

         max_dist = nlist.back().second;

         if (nlist.size() == nfind) {
            if (value > node->GetVarMax() &&
                event.GetDist(node->GetVarMax(), node->GetMod()) > max_dist) {
               return 0;
            }
            if (value < node->GetVarMin() &&
                event.GetDist(node->GetVarMin(), node->GetMod()) > max_dist) {
               return 0;
            }
         }
      }

      const Float_t distance = event.GetDist(node->GetEvent());

      Bool_t insert_this = kFALSE;
      Bool_t remove_back = kFALSE;

      if (nlist.size() < nfind) {
         insert_this = kTRUE;
      }
      else if (nlist.size() == nfind) {
         if (distance < max_dist) {
            insert_this = kTRUE;
            remove_back = kTRUE;
         }
      }
      else {
         std::cerr << "TMVA::kNN::Find() - logic error in recursive procedure" << std::endl;
         return 1;
      }

      if (insert_this) {
         // need typename keyword because qualified dependent names
         // are not valid types unless preceded by 'typename'.
         typename std::list<std::pair<const Node<T> *, Float_t> >::iterator lit = nlist.begin();

         // find a place where current node should be inserted
         for (; lit != nlist.end(); ++lit) {
            if (distance < lit->second) {
               break;
            }
            else {
               continue;
            }
         }

         nlist.insert(lit, std::pair<const Node<T> *, Float_t>(node, distance));

         if (remove_back) {
            nlist.pop_back();
         }
      }
   }

   UInt_t count = 1;
   if (node->GetNodeL() && node->GetNodeR()) {
      if (value < node->GetVarDis()) {
         count += Find(nlist, node->GetNodeL(), event, nfind);
         count += Find(nlist, node->GetNodeR(), event, nfind);
      }
      else {
         count += Find(nlist, node->GetNodeR(), event, nfind);
         count += Find(nlist, node->GetNodeL(), event, nfind);
      }
   }
   else {
      if (node->GetNodeL()) {
         count += Find(nlist, node->GetNodeL(), event, nfind);
      }
      if (node->GetNodeR()) {
         count += Find(nlist, node->GetNodeR(), event, nfind);
      }
   }

   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// This is a global templated function that searches for k-nearest neighbors.
/// list contains all nodes that are closest to event
/// and have sum of event weights >= nfind.
/// Only nodes with positive weights are added to list.
/// Requirement for used classes:
///  - each node contains maximum and minimum values of splitting variable
///    for all its children
///  - min and max range is checked to avoid descending into
///    nodes that are definitely outside current minimum neighbourhood.
///
/// This function should be modified with care.

template<class T>
UInt_t TMVA::kNN::Find(std::list<std::pair<const TMVA::kNN::Node<T> *, Float_t> > &nlist,
                       const TMVA::kNN::Node<T> *node, const T &event, const Double_t nfind, Double_t ncurr)
{

   if (!node || !(nfind < 0.0)) {
      return 0;
   }

   const Float_t value = event.GetVar(node->GetMod());

   if (node->GetWeight() > 0.0) {

      Float_t max_dist = 0.0;

      if (!nlist.empty()) {

         max_dist = nlist.back().second;

         if (!(ncurr < nfind)) {
            if (value > node->GetVarMax() &&
                event.GetDist(node->GetVarMax(), node->GetMod()) > max_dist) {
               return 0;
            }
            if (value < node->GetVarMin() &&
                event.GetDist(node->GetVarMin(), node->GetMod()) > max_dist) {
               return 0;
            }
         }
      }

      const Float_t distance = event.GetDist(node->GetEvent());

      Bool_t insert_this = kFALSE;

      if (ncurr < nfind) {
         insert_this = kTRUE;
      }
      else if (!nlist.empty()) {
         if (distance < max_dist) {
            insert_this = kTRUE;
         }
      }
      else {
         std::cerr << "TMVA::kNN::Find() - logic error in recursive procedure" << std::endl;
         return 1;
      }

      if (insert_this) {
         // (re)compute total current weight when inserting a new node
         ncurr = 0;

         // need typename keyword because qualified dependent names
         // are not valid types unless preceded by 'typename'.
         typename std::list<std::pair<const Node<T> *, Float_t> >::iterator lit = nlist.begin();

         // find a place where current node should be inserted
         for (; lit != nlist.end(); ++lit) {
            if (distance < lit->second) {
               break;
            }

            ncurr += lit -> first -> GetWeight();
         }

         lit = nlist.insert(lit, std::pair<const Node<T> *, Float_t>(node, distance));

         for (; lit != nlist.end(); ++lit) {
            ncurr += lit -> first -> GetWeight();
            if (!(ncurr < nfind)) {
               ++lit;
               break;
            }
         }

         if(lit != nlist.end())
            {
               nlist.erase(lit, nlist.end());
            }
      }
   }

   UInt_t count = 1;
   if (node->GetNodeL() && node->GetNodeR()) {
      if (value < node->GetVarDis()) {
         count += Find(nlist, node->GetNodeL(), event, nfind, ncurr);
         count += Find(nlist, node->GetNodeR(), event, nfind, ncurr);
      }
      else {
         count += Find(nlist, node->GetNodeR(), event, nfind, ncurr);
         count += Find(nlist, node->GetNodeL(), event, nfind, ncurr);
      }
   }
   else {
      if (node->GetNodeL()) {
         count += Find(nlist, node->GetNodeL(), event, nfind, ncurr);
      }
      if (node->GetNodeR()) {
         count += Find(nlist, node->GetNodeR(), event, nfind, ncurr);
      }
   }

   return count;
}

#endif

