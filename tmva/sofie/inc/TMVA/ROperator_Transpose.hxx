#ifndef TMVA_SOFIE_ROPERATOR_TRANSPOSE
#define TMVA_SOFIE_ROPERATOR_TRANSPOSE

#include "TMVA/SOFIE_common.hxx"
#include "TMVA/ROperator.hxx"
#include "TMVA/RModel.hxx"

#include <sstream>

namespace TMVA{
namespace Experimental{
namespace SOFIE{




template <typename T>
class ROperator_Transpose final : public ROperator
{

private:
   std::vector<int_t> fAttrPerm;

   std::string fNData;
   std::string fNOutput;
   std::vector<size_t> fShapeData;
   std::vector<size_t> fShapeOutput;

public:

   ROperator_Transpose() = delete;
   ROperator_Transpose(std::vector<int_t> attr_perm, std::string nameData, std::string nameOutput):
      fAttrPerm(attr_perm), fNData(UTILITY::Clean_name(nameData)), fNOutput(UTILITY::Clean_name(nameOutput)) {
   }

   ROperator_Transpose(std::string nameData, std::string nameOutput):
      fNData(UTILITY::Clean_name(nameData)), fNOutput(UTILITY::Clean_name(nameOutput)) {
   }

   std::vector<ETensorType> TypeInference(std::vector<ETensorType> input){
      return input;
   }

   std::vector<std::vector<size_t>> ShapeInference(std::vector<std::vector<size_t>> input){
      if (input.size() > 1) throw std::runtime_error("TMVA SOFIE Tranpose Op Shape Inference only need 1 input tensor");
      auto& data = input[0];
      std::vector<size_t> output_shape(fAttrPerm.size());
      for (size_t i = 0; i < fAttrPerm.size(); i++){
         output_shape[fAttrPerm[i]] = data[i];
      }
      std::vector<std::vector<size_t>> ret;
      ret.push_back(output_shape);
      return ret;
   }


   void Initialize(RModel& model){
      if (model.CheckIfTensorAlreadyExist(fNData) == false){   //input must be a graph input, or already initialized intermediate tensor
         throw std::runtime_error("TMVA SOFIE Tranpose Op Input Tensor is not found in model");
      }
      fShapeData = model.GetTensorShape(fNData);
      if (fAttrPerm.empty()){
         for (int i = fShapeData.size() - 1; i >= 0; i--){
            fAttrPerm.push_back(i);
         }
      }

      std::vector<size_t> output_shape(fAttrPerm.size());
      for (size_t i = 0; i < fAttrPerm.size(); i++){
         output_shape[fAttrPerm[i]] = fShapeData[i];
      }

      model.AddIntermediateTensor(fNOutput, model.GetTensorType(fNData), output_shape);
      fShapeOutput = output_shape;
   }

   std::string Generate(std::string OpName){
      OpName = "op_" + OpName;
      if (fShapeData.empty() || fShapeOutput.empty()){
         throw std::runtime_error("TMVA SOFIE Transpose Op called to Generate without being initialized first");
      }
      int dim = fShapeData.size();
      int length=1;
      std::vector<int> sizeofindex(dim);
      for (int i = dim - 1; i>=0; i--){
         sizeofindex[i] = length;
         length *= fShapeData[i];
      }
      std::vector<int> index_goto(dim);
      for (int i = 0; i < dim; i++){
         index_goto[fAttrPerm[i]] = i;
      }
      std::vector<int> new_sizeofindex(dim);
      int t = 1;
      for (int i = dim - 1; i>=0; i--){
         new_sizeofindex[i] = t;
         t *= fShapeOutput[i];
      }

      std::stringstream out;
      out << "\t" << "for (int id = 0; id < " << length << " ; id++){\n";
      out << "\t\t " << "tensor_" << fNOutput << "[";
      for (int i =0; i < dim; i++){
         out << "id / " << sizeofindex[i] << " % " << fShapeData[i] << " * " << new_sizeofindex[index_goto[i]];
         if (i != dim - 1) out << " + ";
      }
      out << "] = " << "tensor_" << fNData << "[id];\n";
      out << "\t}\n";
      return out.str();
   }


};

}//SOFIE
}//Experimental
}//TMVA


#endif //TMVA_SOFIE_ROPERATOR_TRANSPOSE
