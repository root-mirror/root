// @(#)root/hist:$Id$
// Author: Olivier Couet 13/07/09

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TVirtualPad.h"
#include "TGraphStruct.h"

#include <cstdio>
#include <iostream>

#include <gvc.h>
#include <gvplugin.h>

#ifdef GVIZ_STATIC
extern gvplugin_library_t gvplugin_dot_layout_LTX_library;
///extern gvplugin_library_t gvplugin_neato_layout_LTX_library;
///extern gvplugin_library_t gvplugin_core_LTX_library;


lt_symlist_t lt_preloaded_symbols[] = {
   { "gvplugin_dot_layout_LTX_library",   (void*)(&gvplugin_dot_layout_LTX_library) },
///   { "gvplugin_neato_layout_LTX_library", (void*)(&gvplugin_neato_layout_LTX_library) },
///   { "gvplugin_core_LTX_library",         (void*)(&gvplugin_core_LTX_library) },
   { 0, 0 }
};
#endif

ClassImp(TGraphStruct);

/** \class TGraphStruct
\ingroup gviz

The Graph Structure is an interface to the graphviz package.

The graphviz package is a graph visualization system. This interface consists in
three classes:

  - TGraphStruct: holds the graph structure. It uses the graphviz library to
    layout the graphs and the ROOT graphics to paint them.
  - TGraphNode: Is a graph node object which can be added in a TGraphStruct.
  - TGraphEdge: Is an edge object connecting two nodes which can be added in
    a TGraphStruct.

Begin_Macro(source)
../../../tutorials/graphs/graphstruct.C
End_Macro

A graph structure can be dumped into a "dot" file using DumpAsDotFile.
*/

////////////////////////////////////////////////////////////////////////////////
/// Graph Structure default constructor.

TGraphStruct::TGraphStruct()
{
   fNodes   = 0;
   fEdges   = 0;
   fGVGraph = 0;
   fGVC     = 0;

   SetMargin();
}

////////////////////////////////////////////////////////////////////////////////
/// Graph Structure default destructor.

TGraphStruct::~TGraphStruct()
{
   gvFreeLayout(fGVC,(Agraph_t*)fGVGraph);
   agclose((Agraph_t*)fGVGraph);
   gvFreeContext(fGVC);

   if (fNodes) delete fNodes;
   if (fEdges) delete fEdges;

}

////////////////////////////////////////////////////////////////////////////////
/// Add the edge "edge" in this TGraphStruct.

void TGraphStruct::AddEdge(TGraphEdge *edge)
{
   if (!fEdges) fEdges = new TList;

   fEdges->Add(edge);
}

////////////////////////////////////////////////////////////////////////////////
/// Create an edge between n1 and n2 and put it in this graph.
///
/// Two edges can connect the same nodes the same way, so there
/// is no need to check if an edge already exists.

TGraphEdge *TGraphStruct::AddEdge(TGraphNode *n1, TGraphNode *n2)
{
   if (!fEdges) fEdges = new TList;

   TGraphEdge *edge = new TGraphEdge(n1, n2);
   fEdges->Add(edge);

   return edge;
}

////////////////////////////////////////////////////////////////////////////////
/// Add the node "node" in this TGraphStruct.

void TGraphStruct::AddNode(TGraphNode *node)
{
   if (!fNodes) fNodes = new TList;

   fNodes->Add(node);
}

////////////////////////////////////////////////////////////////////////////////
/// Create the node "name" if it does not exist and add it to this TGraphStruct.

TGraphNode *TGraphStruct::AddNode(const char *name, const char *title)
{
   if (!fNodes) fNodes = new TList;

   TGraphNode *node = (TGraphNode*)fNodes->FindObject(name);

   if (!node) {
      node = new TGraphNode(name, title);
      fNodes->Add(node);
   }

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Dump this graph structure as a "dot" file.

void TGraphStruct::DumpAsDotFile(const char *filename)
{
   if (!fGVGraph) {
     Int_t ierr = Layout();
     if (ierr) return;
   }
   FILE  *file;
   file=fopen(filename,"wt");
   if (file) {
      agwrite((Agraph_t*)fGVGraph, file);
      fclose(file);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the graph

void TGraphStruct::Draw(Option_t *option)
{
   if (!fGVGraph) {
     Int_t ierr = Layout();
     if (ierr) return;
   }

   // Get the bounding box
   if (gPad) {
      gPad->Range(GD_bb((Agraph_t*)fGVGraph).LL.x-fMargin, GD_bb((Agraph_t*)fGVGraph).LL.y-fMargin,
                  GD_bb((Agraph_t*)fGVGraph).UR.x+fMargin, GD_bb((Agraph_t*)fGVGraph).UR.y+fMargin);
   }

   AppendPad(option);

   // Draw the nodes
   if (fNodes) {
      TGraphNode *node;
      node = (TGraphNode*) fNodes->First();
      node->Draw();
      for(Int_t i = 1; i < fNodes->GetSize(); i++){
         node = (TGraphNode*)fNodes->After(node);
         if (node) node->Draw();
      }
   }

   // Draw the edges
   if (fEdges) {
      TGraphEdge *edge;
      edge = (TGraphEdge*) fEdges->First();
      edge->Draw();
      for(Int_t i = 1; i < fEdges->GetSize(); i++){
         edge = (TGraphEdge*)fEdges->After(edge);
         if (edge) edge->Draw();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Layout the graph into a GraphViz data structure

Int_t TGraphStruct::Layout()
{
   TGraphNode *node;
   TGraphEdge *edge;

   // Create the graph context.
   if (fGVC) gvFreeContext(fGVC);
#ifdef GVIZ_STATIC
   fGVC = gvContextPlugins(lt_preloaded_symbols, 0);
#else
   fGVC = gvContext();
#endif

   // Create the graph.
   if (fGVGraph) {
      gvFreeLayout(fGVC,(Agraph_t*)fGVGraph);
      agclose((Agraph_t*)fGVGraph);
   }
#ifdef WITH_CGRAPH
   fGVGraph = (GVizAgraph_t*)agopen((char*)"GVGraph", Agdirected, 0);
#else
   fGVGraph = (GVizAgraph_t*)agopen((char*)"GVGraph", AGDIGRAPH);
#endif

   // Put the GV nodes into the GV graph
   if (fNodes) {
      node = (TGraphNode*) fNodes->First();
      node->CreateGVNode(fGVGraph);
      for(Int_t i = 1; i < fNodes->GetSize(); i++){
         node = (TGraphNode*)fNodes->After(node);
         if (node) node->CreateGVNode(fGVGraph);
      }
   }

   // Put the edges into the graph
   if (fEdges) {
      edge = (TGraphEdge*) fEdges->First();
      edge->CreateGVEdge(fGVGraph);
      for(Int_t i = 1; i < fEdges->GetSize(); i++){
         edge = (TGraphEdge*)fEdges->After(edge);
         if (edge) edge->CreateGVEdge(fGVGraph);
      }
   }

   // Layout the graph
   int ierr = gvLayout(fGVC, (Agraph_t*)fGVGraph, (char*)"dot");
   if (ierr) return ierr;

   // Layout the nodes
   if (fNodes) {
      node = (TGraphNode*) fNodes->First();
      node->Layout();
      for(Int_t i = 1; i < fNodes->GetSize(); i++){
         node = (TGraphNode*)fNodes->After(node);
         if (node) node->Layout();
      }
   }

   // Layout the edges
   if (fEdges) {
      edge = (TGraphEdge*) fEdges->First();
      edge->Layout();
      for(Int_t i = 1; i < fEdges->GetSize(); i++){
         edge = (TGraphEdge*)fEdges->After(edge);
         if (edge) edge->Layout();
      }
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TGraphStruct::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   out<<"   TGraphStruct *graphstruct = new  TGraphStruct();"<<std::endl;

   // Save the nodes
   if (fNodes) {
      TGraphNode *node;
      node = (TGraphNode*) fNodes->First();
      out<<"   TGraphNode *"<<node->GetName()<<" = graphstruct->AddNode(\""<<
                            node->GetName()<<"\",\""<<
                            node->GetTitle()<<"\");"<<std::endl;
      node->SaveAttributes(out);
      for(Int_t i = 1; i < fNodes->GetSize(); i++){
         node = (TGraphNode*)fNodes->After(node);
         if (node) {
            out<<"   TGraphNode *"<<node->GetName()<<" = graphstruct->AddNode(\""<<
                                  node->GetName()<<"\",\""<<
                                  node->GetTitle()<<"\");"<<std::endl;
            node->SaveAttributes(out);
         }
      }
   }

   // Save the edges
   if (fEdges) {
      TGraphEdge *edge;
      Int_t en = 1;
      edge = (TGraphEdge*) fEdges->First();
      out<<"   TGraphEdge *"<<"e"<<en<<
                            " = new TGraphEdge("<<
                            edge->GetNode1()->GetName()<<","<<
                            edge->GetNode2()->GetName()<<");"<<std::endl;
      out<<"   graphstruct->AddEdge("<<"e"<<en<<");"<<std::endl;
      edge->SaveAttributes(out,Form("e%d",en));
      for(Int_t i = 1; i < fEdges->GetSize(); i++){
         en++;
         edge = (TGraphEdge*)fEdges->After(edge);
         if (edge) {
            out<<"   TGraphEdge *"<<"e"<<en<<
                                  " = new TGraphEdge("<<
                                  edge->GetNode1()->GetName()<<","<<
                                  edge->GetNode2()->GetName()<<");"<<std::endl;
            out<<"   graphstruct->AddEdge("<<"e"<<en<<");"<<std::endl;
            edge->SaveAttributes(out,Form("e%d",en));
         }
      }
   }

   out<<"   graphstruct->Draw();"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////

void TGraphStruct::Streamer(TBuffer &/*b*/)
{
}
