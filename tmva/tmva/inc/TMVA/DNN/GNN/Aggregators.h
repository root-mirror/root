#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/TensorDataLoader.h"
#include "TMVA/DNN/Utility.h"
#include "TMVA/DNN/GNN/GraphData.h"
#include "TMVA/DNN/GNN/Updaters.h"

#ifndef TMVA_DNN_GNN_AGGREGATORS
#define TMVA_DNN_GNN_AGGREGATORS

template<typename Architecture_t>
class GNNEdgeAggregator {

   void Forward(Tensor_t &E, Tensor_t &E_aggregated){
      
      Tensor_t E_aggregated;
      
      bool initialized = false;

      for (int i = 0; i < E.fShape[0]; i++){
         if (initialized) {
            aggregated += E [i];
         }
         
         else {
            e_aggregated = E[i].Copy();
            initialized = true;
         }  
      }

   void Backward (Tensor_t &OriginalE, Tensor_t &GradInput, Tensor_t &GradOutput){
      Tensor_t E_aggregated;
      
      bool initialized = false;

      for (int i = 0; i < E.fShape[0]; i++){
         if (initialized) {
            aggregated += E [i];
         }
         
         else {
            e_aggregated = E[i].Copy();
            initialized = true;
         }  
      }
   }
}


template<typename Architecture_t>
class GNNVerticeAggregator {
   void Forward(Tensor_t &V, Tensor_t &E_aggregated){
      
      Tensor_t V_aggregated;
      
      bool initialized = false;

      for (int i = 0; i < V.fShape[0]; i++){
      if (initialized) {
         V_aggregated += V [i];
      }
      else {
         V_aggregated = V[i].Copy();
         initialized = true;
      }
   
   }
}


template<typename Architecture_t>
class GNNVerticeAggregator {
   void Forward(Tensor_t &V, Tensor_t &E_aggregated){
      
      Tensor_t V_aggregated;
      
      bool initialized = false;

      for (int i = 0; i < V.fShape[0]; i++){
      if (initialized) {
         V_aggregated += V [i];
      }
      else {
         V_aggregated = V[i].Copy();
         initialized = true;
      }
   
   }
}

