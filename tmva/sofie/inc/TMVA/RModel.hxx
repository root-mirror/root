#ifndef TMVA_SOFIE_RMODEL
#define TMVA_SOFIE_RMODEL

#include <vector>
#include <unordered_map>
#include <iostream>
#include <memory>
#include <ctime>
#include <set>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "TMVA/SOFIE_common.hxx"
#include "TMVA/ROperator.hxx"
#include "TBuffer.h"

namespace TMVA{
namespace Experimental{
namespace SOFIE{

class RModel: public TObject{

private:

   std::unordered_map<std::string, InputTensorInfo> fInputTensorInfos; //graph input only; not including operator input (intermediate tensors)
   std::unordered_map<std::string, TensorInfo> fReadyInputTensorInfos;
   std::unordered_map<std::string, InitializedTensor> fInitializedTensors;
   std::unordered_map<std::string, TensorInfo> fIntermediateTensorInfos;
   std::vector<std::string> fOutputTensorNames;

   std::vector<std::unique_ptr<ROperator>> fOperators;

   std::string fName="UnnamedModel";
   std::string fFileName; //file name of original model file for identification
   std::string fParseTime; //UTC date and time string at parsing


   std::string fGC; //generated code
   bool fNeedGemm = true;

   const std::vector<std::string> fAllowedStdLib = {"algorithm"};
   std::set<std::string> fNeededStdLib = {"vector"};



public:

   //explicit move ctor/assn
   RModel(RModel&& other);

   RModel& operator=(RModel&& other);

   //disallow copy
   RModel(const RModel& other) = delete;
   RModel& operator=(const RModel& other) = delete;

   RModel(){}
   RModel(std::string name, std::string parsedtime);

   const std::vector<size_t>& GetTensorShape(std::string name);
   const ETensorType& GetTensorType(std::string name);

   bool CheckIfTensorAlreadyExist(std::string tensor_name);
   void AddInputTensorInfo(std::string input_name, ETensorType type, std::vector<Dim> shape);
   void AddInputTensorInfo(std::string input_name, ETensorType type, std::vector<size_t> shape);
   void AddOperator(std::unique_ptr<ROperator> op, int order_execution = -1);
   void AddInitializedTensor(std::string tensor_name, ETensorType type, std::vector<std::size_t> shape, std::shared_ptr<void> data);
   void AddIntermediateTensor(std::string tensor_name, ETensorType type, std::vector<std::size_t> shape);
   void AddNeededStdLib(std::string libname){
      for (auto& i: fAllowedStdLib){
         if ( i == libname) fNeededStdLib.insert(libname);
      }
   }
   void AddOutputTensorNameList(std::vector<std::string> outputtensornames);
   void UpdateInitializedTensor(std::string tensor_name, ETensorType type, std::vector<std::size_t> shape, std::shared_ptr<void> data);
   std::shared_ptr<void> GetInitializedTensorData(std::string tensor_name);


   void Initialize();
   void Generate();

   void PrintGenerated(){
      std::cout << fGC;
   }
   void PrintIntermediateTensors();
   void OutputGenerated(std::string filename = "");


/*
   template <typename T>
   void AddInitializedTensor(std::string tensor_name, RTensor<T> new_tensor){
      //a view only
      T obj;
      if (fInitializedTensors.find(tensor_name) != fInitializedTensors.end()){
         throw std::runtime_error("TMVA-SOFIE: initialized tensor with name " + tensor_name + " already exists \n");
      }
      InitializedTensor new_tensor_ {GetTemplatedType(obj), new_tensor.GetShape() , static_cast<void>(new_tensor.GetData())};
      fInitializedTensors[tensor_name] = new_tensor_;
   }
*/

   void PrintRequiredInputTensors();
   void PrintInitializedTensors();
   void HeadInitializedTensors(std::string name, int n_print = 50);

   ~RModel(){
      /*
      for (auto& i: fInitializedTensors){
         free(i.second.data);
      }
      */
   }
   ClassDef(RModel,1);
};

}//SOFIE
}//Experimental
}//TMVA

#endif //TMVA_SOFIE_RMODEL
