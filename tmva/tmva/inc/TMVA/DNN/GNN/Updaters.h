#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/TensorDataLoader.h"
#include "TMVA/DNN/Utility.h"
#include "TMVA/DNN/GNN/GraphData.h"
#include "TMVA/DNN/GNN/Aggregators.h"


#ifndef TMVA_DNN_GNN_UPDATERS
#define TMVA_DNN_GNN_UPDATERS


template<typename Architecture_t>
class GNNEdgeUpdater{
public:

   using Tensor_t = typename Architecture_t::Tensor_t;
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;
   using EdgeList = std::vector <std::pair <std::size, std::size>>;
   using Shape_t = std::vector<std::size_t>;

private:
   TDeepNet <Architecture_t> updater_network;
   EOutputFunction f;
   Shape_t flatten_shape;
   size_t E_size;
   size_t V_size;
   size_t G_size;
   Shape_t E_shape;
   Shape_t V_shape;
   Shape_t G_shape;
   
public:

   GNNEdgeUpdater(TDeepNet <Architecture_t> updater_network, EOutputFunction f){
      this->updater_network = updater_network;
      this->f = f;
   }

   void Forward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t & output, EdgeList &E_connectivity_updated){
      nEdges = E.GetNrows();
      Tensor_t new_e(E.fShape);
      for (int i = 0; i < nEdges; i++){

         std::size source_v_id = E_connectivity[i].first;
         std::size targer_v_id = E_connectivity[i].second;

         Tensor_t flatten_v_s = V[source_v_id].Flatten();
         Tensor_t flatten_v_t = V[targer_v_id].Flatten();
         
         Tensor_t flatten_e = E[i].Flatten();
         Tensor_t flatten_g = G.Flatten();

         Tensor_t x = flatten_v_s.Concatenate(flatten_v_t).Concatenate(flatten_e).Concatenate(flatten_g);

         flatten_shape = x.fShape();
      
         updater_network.Forward(x)

         Matrix_t preditcion_matrix (E[i].fSize, 1);
         updater_network.Prediction(preditcion_matrix, f);

         Tensor_t res (preditcion_matrix);
         
         new_e[i].CopyData(res);
      }
      this -> E_size = E.fSize();
      this -> V_size = V.fSize();
      this -> Z_size = Z.fSize();
      this -> E_shape = E.fShape();
      this -> V_shape = V.fShape();
      this -> Z_shape = Z.fShape();

      E_connectivity_updated = E_connectivity
      output = new_e;
   }
/*
   void Backward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t & correct_output, Matrix_t &E_grad, Matrix_t &V_grad, Matrix_t &G_grad){
      Matrix_t weights(1, 1);
      fillMatrix(weights, 1.0);
      
      nEdges = E.GetNrows();
      Tensor_t new_e(flatten_shape);
      
      Tensor_t grad_t(flatten_shape);
      Matrix_t global_grad = grad_t.GetMatrix();
      fillMatrix(global_grad, 0.0);
      
      Tensor_t grad_div_t(flatten_shape);
      Matrix_t grad_div = grad_div_t.GetMatrix();
      fillMatrix(grad_div, nEdges);
      
      for (int i = 0; i < nEdges; i++){
         std::size source_v_id = E_connectivity[i].first;
         std::size targer_v_id = E_connectivity[i].second;

         Tensor_t flatten_v_s = V[source_v_id].Flatten();
         Tensor_t flatten_v_t = V[targer_v_id].Flatten();
         
         Tensor_t flatten_e = E[i].Flatten();
         Tensor_t flatten_g = G.Flatten();

         Tensor_t x = flatten_v_s.Concatenate(flatten_v_t).Concatenate(flatten_e).Concatenate(flatten_g);

         Matrix_t grad;

         updater_network.Backward(x, correct_output.GetMatrix(), weights, grad);

         global_grad += grad;
      }

      global_grad /= grad_div;

      Tensor_t global_grad_tensor (global_grad);
      data.

   }
*/
}



template<typename Architecture_t>
class GNNNodeUpdater{
public:

   using Tensor_t = typename Architecture_t::Tensor_t;
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;
   using EdgeList = std::vector <std::pair <std::size, std::size>>;
   
private:
   TDeepNet <Architecture_t> updater_network;
   EOutputFunction f;

   
public:

   GNNEdgeUpdater(TDeepNet <Architecture_t> updater_network, EOutputFunction f){
      this->updater_network = updater_network;
      this->f = f;
   }

   void Forward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t & output, GNNEdgePerVerticeAggregator <Architecture_t> &edge_per_vertice_aggregator){
      Tensor_t new_v(V.fShape);
    

      nVertices = V.GetNrows();
      Tensor_t new_v(V.fShape);
      for (int i = 0; i < nVertices; i++){

         Tensor_t e_aggregated = edge_per_vertice_aggregator(Tensor_t E, EdgeList E_connectivity);

         Tensor_t flatten_v = vertices[i].Flatten();      
         Tensor_t flatten_e = e_aggregated.Flatten();
         Tensor_t flatten_g = global.Flatten();

         Tensor_t x = flatten_v.Concatenate(flatten_e).Concatenate(flatten_g);

         updater_network.Forward(x)

         Matrix_t preditcion_matrix (V[i].fSize, 1);
         updater_network.Prediction(preditcion_matrix, f);

         Tensor_t res (preditcion_matrix);
         
         new_v[i].CopyData(res);
      }

      output = new_v;

   }

/*
   void Backward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t & correct_output, Matrix_t &activationGradients){
      Matrix_t weights(1, 1);
      fillMatrix(weights, 1.0);
      
      nVertices = V.GetNrows();
      Tensor_t new_e(V.fShape);
      
      Tensor_t grad_t(V.fShape);
      Matrix_t global_grad = grad_t.GetMatrix();
      fillMatrix(global_grad, 0.0);
      
      Tensor_t grad_div_t(V.fShape);
      Matrix_t grad_div = grad_div_t.GetMatrix();
      fillMatrix(grad_div, nEdges);
      
      for (int i = 0; i < nEdges; i++){

         std::size source_v_id = E_connectivity[i].first;
         std::size targer_v_id = E_connectivity[i].second;

         Tensor_t flatten_v_s = V[source_v_id].Flatten();
         Tensor_t flatten_v_t = V[targer_v_id].Flatten();
         
         Tensor_t flatten_e = E[i].Flatten();
         Tensor_t flatten_g = G.Flatten();

         Tensor_t x = flatten_v_s.Concatenate(flatten_v_t).Concatenate(flatten_e).Concatenate(flatten_g);

         Matrix_t grad;

         updater_network.Backward(x, correct_output.GetMatrix(), weights, grad);

         global_grad += grad;
      }

      global_grad /= grad_div;
      activationGradients = global_grad;

   }
*/
}



template<typename Architecture_t>
class GNNGlobalUpdater{
public:

   using Tensor_t = typename Architecture_t::Tensor_t;
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;
   using EdgeList = std::vector <std::pair <std::size, std::size>>;
   
private:
   TDeepNet <Architecture_t> updater_network;
   EOutputFunction f;

   
public:

   GNNEdgeUpdater(TDeepNet <Architecture_t> updater_network, EOutputFunction f){
      this->updater_network = updater_network;
      this->f = f;
   }

   void Forward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t &new_G, GNNEdgeAggregator <Architecture_t> &edge_aggregator, GNNVerticeAggregator <Architecture_t> &vertice_aggregator){
      Tensor_t new_g(G.fShape);

      Tensor_t e_aggregated = edge_aggregator(Tensor_t E);
      Tensor_t v_aggregated = vertice_aggregator(Tensor_t V);

      Tensor_t flatten_v = v_aggregated.Flatten();      
      Tensor_t flatten_e = e_aggregated.Flatten();
      Tensor_t flatten_g = global.Flatten();

      Tensor_t x = flatten_v.Concatenate(flatten_e).Concatenate(flatten_g);

      updater_network.Forward(x)

      Matrix_t preditcion_matrix (G[i].fSize, 1);
      updater_network.Prediction(preditcion_matrix, f);

      new_G = Tensor_t res (preditcion_matrix);
            
   }

/*
   void Backward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t & correct_output, Matrix_t &activationGradients){
      

   }
*/
}





template<typename Architecture_t>
class GNNEdgeUpdaterIdentity : GNNEdgeUpdater{
public:

   GNNEdgeUpdaterIdentity(){};

   void Forward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t & output){
      return output = E;
   }

   void Backward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t & correct_output, Matrix_t &activationGradients){
      Tensor_t grad_t(E.fShape);
      Matrix_t global_grad = grad_t.GetMatrix();
      fillMatrix(global_grad, 0.0);
      activationGradients = global_grad;
   }
}



template<typename Architecture_t>
class GNNVerticeUpdaterIdentity : GNNVerticeUpdater{
public:

   GNNVerticeUpdaterIdentity(){};

   void Forward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t & output){
      return output = V;
   }

   void Backward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t & correct_output, Matrix_t &activationGradients){
      Tensor_t grad_t(V.fShape);
      Matrix_t global_grad = grad_t.GetMatrix();
      fillMatrix(global_grad, 0.0);
      activationGradients = global_grad;
   }
}

template<typename Architecture_t>
class GNNGlobalUpdaterIdentity : GNNGlobalUpdater{
public:

   GNNGlobalUpdaterIdentity(){};

   void Forward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t & output){
      return output = G;
   }

   void Backward(Tensor_t &E, Tensor_t &V, Tensor_t &G, EdgeList &E_connectivity, Tensor_t & correct_output, Matrix_t &activationGradients){
      Tensor_t grad_t(G.fShape);
      Matrix_t global_grad = grad_t.GetMatrix();
      fillMatrix(global_grad, 0.0);
      activationGradients = global_grad;
   }
}



