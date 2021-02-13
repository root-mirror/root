// Author: Enrico Guiraud, Danilo Piparo, CERN, Massimo Tumolo Politecnico di Torino 08/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_GRAPHUTILS
#define ROOT_GRAPHUTILS

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <type_traits>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDF/GraphNode.hxx>

#include <iostream>

namespace ROOT {
namespace Detail {
namespace RDF {
class RDefineBase;
class RFilterBase;
class RRangeBase;
} // namespace RDF
} // namespace Detail

namespace Internal {
namespace RDF {
class RBookedDefines;

namespace GraphDrawing {
std::shared_ptr<GraphNode>
CreateDefineNode(const std::string &columnName, const ROOT::Detail::RDF::RDefineBase *columnPtr);

std::shared_ptr<GraphNode> CreateFilterNode(const ROOT::Detail::RDF::RFilterBase *filterPtr);

std::shared_ptr<GraphNode> CreateRangeNode(const ROOT::Detail::RDF::RRangeBase *rangePtr);


/// Add the Defines that have been added between this node and the previous to the graph.
/// Return the new "upmost" node, i.e. the last of the Defines added if any, otherwise the node itself
std::shared_ptr<GraphNode> AddDefinesToGraph(std::shared_ptr<GraphNode> node,
                                             const RDFInternal::RBookedDefines &defines,
                                             const std::vector<std::string> &prevNodeDefines);

// clang-format off
/**
\class ROOT::Internal::RDF::GraphCreatorHelper
\ingroup dataframe
\brief Helper class that provides the operation graph nodes.

 This class is the single point from which graph nodes can be retrieved. Every time an object is created,
 it clears the static members and starts again.
 By asking this class to create a node, it will return an existing node if already created, otherwise a new one.
*/
// clang-format on
class GraphCreatorHelper {
private:
   using DefinesNodesMap_t = std::map<const ROOT::Detail::RDF::RDefineBase *, std::weak_ptr<GraphNode>>;
   using FiltersNodesMap_t = std::map<const ROOT::Detail::RDF::RFilterBase *, std::weak_ptr<GraphNode>>;
   using RangesNodesMap_t = std::map<const ROOT::Detail::RDF::RRangeBase *, std::weak_ptr<GraphNode>>;

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Stores the columns defined and which node in the graph defined them.
   static DefinesNodesMap_t &GetStaticColumnsMap()
   {
      static DefinesNodesMap_t sMap;
      return sMap;
   };

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Stores the filters defined and which node in the graph defined them.
   static FiltersNodesMap_t &GetStaticFiltersMap()
   {
      static FiltersNodesMap_t sMap;
      return sMap;
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Stores the ranges defined and which node in the graph defined them.
   static RangesNodesMap_t &GetStaticRangesMap()
   {
      static RangesNodesMap_t sMap;
      return sMap;
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Invoked by the RNodes to create a define graph node.
   friend std::shared_ptr<GraphNode>
   CreateDefineNode(const std::string &columnName, const ROOT::Detail::RDF::RDefineBase *columnPtr);

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Invoked by the RNodes to create a Filter graph node.
   friend std::shared_ptr<GraphNode> CreateFilterNode(const ROOT::Detail::RDF::RFilterBase *filterPtr);

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Invoked by the RNodes to create a Range graph node.
   friend std::shared_ptr<GraphNode> CreateRangeNode(const ROOT::Detail::RDF::RRangeBase *rangePtr);

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Starting from any leaf (Action, Filter, Range) it draws the dot representation of the branch.
   std::string FromGraphLeafToDot(std::shared_ptr<GraphNode> leaf);

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Starting by an array of leaves, it draws the entire graph.
   std::string FromGraphActionsToDot(std::vector<std::shared_ptr<GraphNode>> leaves);

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Starting from the root node, prints the entire graph.
   std::string RepresentGraph(ROOT::RDataFrame &rDataFrame);

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Starting from the root node, prints the entire graph.
   std::string RepresentGraph(RLoopManager *rLoopManager);

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Starting from a Filter or Range, prints the branch it belongs to
   template <typename Proxied, typename DataSource>
   std::string RepresentGraph(RInterface<Proxied, DataSource> &rInterface)
   {
      auto loopManager = rInterface.GetLoopManager();
      loopManager->Jit();

      return FromGraphLeafToDot(rInterface.GetProxiedPtr()->GetGraph());
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Starting from an action, prints the branch it belongs to
   template <typename T>
   std::string RepresentGraph(const RResultPtr<T> &resultPtr)
   {
      auto loopManager = resultPtr.fLoopManager;
      if (!loopManager)
         throw std::runtime_error("Something went wrong");

      if (std::is_same<T, RInterface<RLoopManager, void>>::value) {
         return RepresentGraph(loopManager);
      }

      loopManager->Jit();

      auto actionPtr = resultPtr.fActionPtr;
      return FromGraphLeafToDot(actionPtr->GetGraph());
   }

public:
   ////////////////////////////////////////////////////////////////////////////
   /// \brief Functor. Initializes the static members and delegates the work to the right override.
   /// \tparam NodeType the RNode from which the graph has to be drawn
   template <typename NodeType>
   std::string operator()(NodeType &node)
   {
      // First all static data structures are cleaned, to avoid undefined behaviours if more than one Represent is
      // called
      GetStaticFiltersMap() = FiltersNodesMap_t();
      GetStaticColumnsMap() = DefinesNodesMap_t();
      GetStaticRangesMap() = RangesNodesMap_t();
      GraphNode::ClearCounter();
      // The Represent can now start on a clean environment
      return RepresentGraph(node);
   }
};

} // namespace GraphDrawing
} // namespace RDF
} // namespace Internal
} // namespace ROOT

#endif
