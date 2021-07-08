#include "TMVA/DNN/Architectures/Cpu/CpuBuffer.h"
#include "TMVA/DNN/Architectures/Cpu/CpuMatrix.h"
#include "TMVA/DNN/Architectures/Cpu/CpuTensor.h"
#include "TMVA/DNN/GNN/GraphData.h"
#include "TMVA/DNN/GNN/Updaters.h"


//////////////////////////////////////////////////////////////////////
// This class implements the general structure of a GNN block 
// according the framework proposed in the paper by Battaglia et al.
// "Relational inductive biases, deep learning, and graph networks" by
// P Battaglia et al. (2018)
//////////////////////////////////////////////////////////////////////

#ifndef TMVA_DNN_GNN_BLOCK
#define TMVA_DNN_GNN_BLOCK


namespace TMVA{
namespace DNN{
namespace GNN{

template <typename Architecture_t>
class TGNNBlock : public VGeneralLayer<Architecture_t> {
public:

   using Tensor_t = typename Architecture_t::Tensor_t;
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;
   using EdgeList = std::vector <std::pair <std::size, std::size>>;

private: 
   GNNEdgeUpdater <Architecture_t>  edge_updater;
   GNNVerticeUpdater <Architecture_t>  vertice_updater;
   GNNGlobalUpdater <Architecture_t>  global_updater;
   GNNEdgeAggregator <Architecture_t> edge_aggregator;
   GNNEdgePerVerticeAggregator <Architecture_t> edge_per_vertice_aggregator;
   GNNVerticeAggregator <Architecture_t> vertice_aggregator;
   TGraphData <Architecture_t> last_updated_graph;

public:

   // General constructor
   TGNNBlock  (GNNEdgeUpdater <Architecture_t>  edge_updater,
               GNNVerticeUpdater <Architecture_t>  vertice_updater,
               GNNGlobalUpdater <Architecture_t>  global_updater,
               GNNEdgeAggregator <Architecture_t> edge_aggregator,
               GNNEdgePerVerticeAggregator <Architecture_t> edge_per_vertice_aggregator,
               GNNVerticeAggregator <Architecture_t> vertice_aggregator){
                  this -> edge_updater = edge_updater;
                  this -> vertice_updater = vertice_updater;
                  this -> global_updater = global_updater;
                  this -> edge_aggregator = edge_aggregator;
                  this -> edge_per_vertice_aggregator = edge_per_vertice_aggregator;
                  this -> vertice_aggregator = vertice_aggregator;
               }
   
   // A method implementing the forward pass through the whole block
   void Forward (TGraphData<Architecture_t> &X);
};  


template<typename Architecture_t>
void TGNNBlock<Architecture_t>::Forward (TGraphData <Architecture_t> &X, TGraphData <Architecture_t> &X_updated){
   EdgeList E_connectivity = X.GetEdgeConnectivity();
   Tensor_t V = X.GetVerticeAttributes();
   Tensor_t E = X.GetEdgeAttributes();
   Tensor_t G = X.GetGlobalAttributes();

   EdgeList E_connectivity_updated;
   Tensor_t E_updated;
   Tensor_t V_updated;
   Tensor_t G_updated;

   edge_updater.Forward (E, V, G, E_connectivity, E_updated, E_connectivity_updated);
   vertice_updater.Forward (E_updated, V, G, E_connectivity_updated, V_updated, edge_per_vertice_aggregator);
   global_updater.Forward(E_updated, V_updated, G, E_connectivity_updated, G_updated, edge_per_vertice_aggregator);

   TGraphData <Architecture_t> x_new (E_updated, V_updated, G_updated, E_connectivity_updated);
   X_updated = x_new;
   this -> last_updated_graph = X_updated;
}

void TGNNBlock<Architecture_t>::Backward (TGraphData <Architecture_t> &X, TGraphData <Architecture_t> &X_updated){
   EdgeList E_connectivity = X.GetEdgeConnectivity();
   Tensor_t V = X.GetVerticeAttributes();
   Tensor_t E = X.GetEdgeAttributes();
   Tensor_t G = X.GetGlobalAttributes();

   EdgeList E_connectivity_updated;
   Tensor_t E_updated;
   Tensor_t V_updated;
   Tensor_t G_updated;

   global_updater.Backward(E_updated, V_updated, G, E_connectivity_updated, G_updated, edge_per_vertice_aggregator){
      /*
      I've deleted backward pass, since I believe it was incorrectly implemented
      */
   }

}


} // namespace RNN
} // namespace DNN
} // namespace GNN

#endif

/*
old forward function 

TODO: Clean this, disassembe for parts until August

Tensor_t edges = X.GetEdgeAttributes();
   Tensor_t edgeConnectivity = X.GetEdgeConnectivity();
   EdgeList vertices = X.GetVerticeAttributes();
   Tensor_t global = X.GetGlobalAttributes();

   nEdges = Edges.GetNrows();
   nVertices = Vertices.GetNrows();

   TGraphData <Architecture_t> newX = TGraphData <Architecture_t>(X);

   // Edge attribute update step
   
   Tensor_t new_e(Edges.fShape);
   for (int i = 0; i < nEdges; i++){

      std::size source_v_id = edgeConnectivity[i].first;
      std::size targer_v_id = edgeConnectivity[i].second;

      Tensor_t flatten_v_s = vertices[source_v_id].Flatten();
      Tensor_t flatten_v_t = vertices[targer_v_id].Flatten();
      
      Tensor_t flatten_e = edges[i].Flatten();

      Tensor_t flatten_g = global.Flatten();

      Tensor_t x = flatten_v_s.Concatenate(flatten_v_t).Concatenate(flatten_e).Concatenate(flatten_g);

      new_e[i].CopyData(x);
   }

   // Vertice update 
   Tensor_t new_v(vertices.fShape);
   for (int i = 0; i < nVertices; i++){

      Tensor_t e_aggregated;
      bool e_initialized = false;

      for (int i = 0; i < nEdges; i++){
         if (edgeConnectivity[i].second == i){
            if (e_initialized) {
               e_aggregated += new_e [i];
            }
            else {
               e_aggregated new_e[i].Copy();
               e_initialized;
            }
         }
      }

      Tensor_t flatten_v = vertices[i].Flatten();      
      Tensor_t flatten_e_ag = e_aggregated.Flatten();
      Tensor_t flatten_g = global.Flatten();

      Tensor_t x = flatten_v_s.Concatenate(flatten_v_t).Concatenate(flatten_e).Concatenate(flatten_g);

      vertices[i].CopyData(x);
   }

   // Global update
   Tensor_t new_g(global.fShape);
   
   Tensor_t e_aggregated;
   bool e_initialized = false;
   for (int i = 0; i < nEdges; i++){
      if (e_initialized) {
         e_aggregated += new_e [i];
      }
      else {
         e_aggregated new_e[i].Copy();
         e_initialized;
      }
      
   }

   Tensor_t v_aggregated;
   bool v_initialized = false;
   for (int i = 0; i < nEdges; i++){
      if (v_initialized) {
         v_aggregated += new_v [i];
      }
      else {
         v_aggregated new_v[i].Copy();
         v_initialized;
      }
   
   }

   Tensor_t x = flatten_v_s.Concatenate(global.Flatten()).Concatenate(e_aggregated.Flatten()).Concatenate(v_aggregated.Flatten());


*/