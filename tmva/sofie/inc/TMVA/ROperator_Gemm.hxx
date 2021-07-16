#ifndef TMVA_SOFIE_ROPERATOR_GEMM
#define TMVA_SOFIE_ROPERATOR_GEMM


#include "TMVA/SOFIE_common.hxx"
#include "TMVA/ROperator.hxx"
#include "TMVA/RModel.hxx"

#include <sstream>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <limits>

namespace TMVA{
namespace Experimental{
namespace SOFIE{


   template <typename T>
   class ROperator_Gemm final : public ROperator
   {

   private:
      float fAttrAlpha = 1.0;
      float fAttrBeta = 1.0;
      int_t fAttrTransA = 0;
      int_t fAttrTransB = 0;

      std::string fNA;
      std::string fNB;
      std::string fNC = "";
      std::string fNY;
      std::vector<size_t> fShapeA;
      std::vector<size_t> fShapeB;
      std::vector<size_t> fShapeC;
      std::vector<size_t> fShapeY;

      std::string fType;

   public:

      ROperator_Gemm() = delete;
      ROperator_Gemm(float alpha, float beta, int_t transA, int_t transB, std::string nameA, std::string nameB, std::string nameY):
         fAttrAlpha(alpha), fAttrBeta(beta), fAttrTransA(transA), fAttrTransB(transB), fNA(UTILITY::Clean_name(nameA)),
         fNB(UTILITY::Clean_name(nameB)), fNY(UTILITY::Clean_name(nameY)) {

         if (std::is_same<T, float>::value) {
            fType = "float";
         }else{
            throw std::runtime_error("TMVA SOFIE Encountered unsupported type parsing a gemm operator");
         }
      }

      ROperator_Gemm(float alpha, float beta, int_t transA, int_t transB, std::string nameA, std::string nameB, std::string nameC, std::string nameY):
         fAttrAlpha(alpha), fAttrBeta(beta), fAttrTransA(transA), fAttrTransB(transB), fNA(UTILITY::Clean_name(nameA)),
         fNB(UTILITY::Clean_name(nameB)), fNC(UTILITY::Clean_name(nameC)), fNY(UTILITY::Clean_name(nameY)) {

         if (std::is_same<T, float>::value) {
            fType = "float";
         }else{
            throw std::runtime_error("TMVA SOFIE Encountered unsupported type parsing a gemm operator");
         }
      }

      std::vector<ETensorType> TypeInference(std::vector<ETensorType> input){
         ETensorType out = input[0];
         return {out};
      }

      std::vector<std::vector<size_t>> ShapeInference(std::vector<std::vector<size_t>> input){
         if (input.size() > 3) throw std::runtime_error("TMVA SOFIE Gemm Op Shape Inference only need 2 or 3 input tensor");
         for (auto& i: input){
            if (i.size() > 2){
               throw std::runtime_error("TMVA SOFIE Gemm Op Shape Inference only accept input tensor with 2 dimensions");
            }
         }
         std::vector<std::vector<size_t>> ret;
         if (input.size() == 3){
            ret.push_back(input[2]);   //shape of C is shape of Y
            return ret;
         }
         std::vector<size_t> s_a(input[0]);
         std::vector<size_t> s_b(input[1]);
         if (fAttrTransA){
            std::reverse(s_a.begin(), s_a.end());
         }
         if (fAttrTransB){
            std::reverse(s_b.begin(), s_b.end());
         }
         std::vector<size_t> s_y(2);
         s_y[0] = s_a[0];
         s_y[1] = s_b[1];
         ret.push_back(s_y);
         return ret;
      }



      void Initialize(RModel& model){
         //TODO: propagate A or B as specified by ONNX standard

         if ((model.CheckIfTensorAlreadyExist(fNA) == false) || (model.CheckIfTensorAlreadyExist(fNB) == false) ){   //input must be a graph input, or already initialized intermediate tensor
            throw std::runtime_error("TMVA SOFIE Gemm Op Input Tensor " + fNA + " or " + fNB + " is not found in model");
         }
         if (fNC != ""){
            if (model.CheckIfTensorAlreadyExist(fNC) == false){   //input must be a graph input, or already initialized intermediate tensor
               throw std::runtime_error("TMVA SOFIE Gemm Op Input Tensor" + fNC + " is not found in model");
            }
         }
         fShapeA = model.GetTensorShape(fNA);
         if (fShapeA.size() != 2){
            throw std::runtime_error("TMVA SOFIE Gemm Op Input Tensor" + fNA + " is not of 2 dimensions");
         }
         fShapeB = model.GetTensorShape(fNB);
         if (fShapeB.size() != 2){
            throw std::runtime_error("TMVA SOFIE Gemm Op Input Tensor" + fNB + " is not of 2 dimensions");
         }
         fShapeY = ShapeInference({fShapeA, fShapeB})[0];
         if (fNC != ""){
            fShapeC = model.GetTensorShape(fNC);

            bool broadcast_needed = false;
            for (size_t i =0; i < fShapeC.size(); i++){
               if (fShapeC[i]!=fShapeY[i]){
                  broadcast_needed = true;
                  break;
               }
            }

            if (broadcast_needed){
               auto original_data = model.GetInitializedTensorData(fNC);
               if (fType == "float"){

                  std::shared_ptr<void> new_data_ptr(UTILITY::Unidirectional_broadcast<float>(static_cast<float*>(original_data.get()), fShapeC, fShapeY), std::default_delete<float[]>());


                  model.UpdateInitializedTensor(fNC, model.GetTensorType(fNC), fShapeY, new_data_ptr);
                  fShapeC = fShapeY;
               }
            }
         }




         model.AddIntermediateTensor(fNY, model.GetTensorType(fNA), fShapeY);
         model.AddNeededStdLib("algorithm");

      }



      std::string Generate(std::string OpName){
         OpName = "op_" + OpName;
         if (fShapeA.empty() || fShapeB.empty() || fShapeY.empty() || (fNC != "" && fShapeC.empty())){
            throw std::runtime_error("TMVA SOFIE Gemm Op called to Generate without being initialized first");
         }
         std::stringstream out;
         out <<"\t" << "char " << OpName << "_transA = " << (fAttrTransA ? "\'t\'" : "\'n\'") << ";\n";
         out <<"\t" << "char " << OpName << "_transB = " << (fAttrTransB ? "\'t\'" : "\'n\'") << ";\n";
         int m = (fAttrTransA ? fShapeA[1] : fShapeA[0]);
         int n = (fAttrTransB ? fShapeB[0] : fShapeB[1]);
         int k = (fAttrTransA ? fShapeA[0] : fShapeA[1]);
         out <<"\t" << "int " << OpName << "_m = " << m << ";\n";
         out <<"\t" << "int " << OpName << "_n = " << n << ";\n";
         out <<"\t" << "int " << OpName << "_k = " << k << ";\n";
         out <<"\t" << "float " << OpName << "_alpha = " << std::setprecision(std::numeric_limits<float>::max_digits10) << fAttrAlpha << ";\n";
         out <<"\t" << "float " << OpName << "_beta = " << std::setprecision(std::numeric_limits<float>::max_digits10) << fAttrBeta << ";\n";
         out <<"\t" << "int " << OpName << "_lda = " << (fAttrTransA ? m : k) << ";\n";
         out <<"\t" << "int " << OpName << "_ldb = " << (fAttrTransB ? k : n) << ";\n";
         if (fNC != ""){
            int length = 1;
            for (auto& i: fShapeC){
               length *= i;
            }
            out << "\t" << "std::copy(" << "tensor_" << fNC << ", " << "tensor_" << fNC << " + " << length << ", " << "tensor_" << fNY << ");\n";
         }
         if (fType == "float"){
            out << "\t" << "BLAS::sgemm_(&" << OpName << "_transB, &" << OpName << "_transA, &" << OpName
             << "_n, &" << OpName << "_m, &" << OpName << "_k, &" << OpName << "_alpha, " << "tensor_" << fNB
             << ", &" << OpName << "_ldb, " << "tensor_" << fNA << ", &" << OpName << "_lda, &" << OpName << "_beta, " << "tensor_" << fNY << ", &"
             << OpName << "_n);\n";
          }

          return out.str();

         }



   };


}//SOFIE
}//Experimental
}//TMVA


#endif //TMVA_SOFIE_ROPERATOR_GEMM
