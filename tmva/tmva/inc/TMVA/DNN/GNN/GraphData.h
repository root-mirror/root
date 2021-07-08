#include "TMVA/DNN/Architectures/Cpu/CpuBuffer.h"
#include "TMVA/DNN/Architectures/Cpu/CpuMatrix.h"
#include "TMVA/DNN/Architectures/Cpu/CpuTensor.h"

#ifndef TMVA_DNN_GNN_GRAPHDATA
#define TMVA_DNN_GNN_GRAPHDATA

namespace TMVA{
namespace DNN{
namespace GNN{

template<typename AReal = Float_t>
class GraphData{
   using Tensor_t = TCpuTensor<AReal>;
   using EdgeList = std::vector <std::pair <std::size, std::size>>;

   private:
   // Edges attribute tensor of a shape [N_e, ...]
   Tensor_t EdgeAttributes;

   // Vertices attribute tensor of a shape [N_v, ...]
   Tensor_t VerticeAttributes;

   // Global attribute tensor
   Tensor_t GlobalAttributes;

   // Vector of pairs of edge indices
   EdgeList EdgeConnectivity;

   public:
   GraphData(Tensor_t EdgeAttributes, Tensor_t VerticeAttributes, Tensor_t GlobalAttributes, EdgeList EdgeConnectivity){
      this -> EdgeAttributes = EdgeAttributes;
      this -> VerticeAttributes = VerticeAttributes;
      this -> GlobalAttributes = GlobalAttributes;
      this -> EdgeConnectivity = EdgeConnectivity;
   }

   void SetEdgeAttributes(Tensor_t & t){
      this->EdgeAttributes = t;
   } 

   void SetVerticeAttributes(Tensor_t & t){
      this->VerticeAttributes = t;
   }

   void SetGlobalAttributes(Tensor_t & t){
      this->GlobalAttributes = t;
   }

   void SetEdgeConnectivity(EdgeList & e){
      this->EdgeConnectivity = e;
   }


   Tensor_t & GetEdgeAttributes(){
      return this->EdgeAttributes;
   } 

   Tensor_t & GetVerticeAttributes(){
      return this->VerticeAttributes;
   }

   Tensor_t & GetGlobalAttributes(){
      return this->GlobalAttributes ;
   }

   EdgeList & GetEdgeConnectivity(){
      return this->EdgeConnectivity;
   }

}

}
}
}