#include "TMVA/RModelParser_ONNX.hxx"
#include "onnx_proto3.pb.h"

#include <string>
#include <memory>

namespace TMVA{
namespace Experimental{
namespace SOFIE{

namespace INTERNAL{

std::unique_ptr<ROperator> make_ROperator(size_t idx, const onnx::GraphProto& graphproto, std::unordered_map<std::string, ETensorType>& tensor_type){
   const auto& nodeproto = graphproto.node(idx);
   auto find = mapOptypeOperator.find(nodeproto.op_type());
   if (find == mapOptypeOperator.end()){
      throw std::runtime_error("TMVA::SOFIE - Operator type " + nodeproto.op_type() + " is not yet supported");
   }else{
      return (find->second)(nodeproto, graphproto, tensor_type);
   }
}

std::unique_ptr<ROperator> make_ROperator_Transpose(const onnx::NodeProto& nodeproto, const onnx::GraphProto& graphproto, std::unordered_map<std::string, ETensorType>& tensor_type){

   ETensorType input_type;

   auto input_name =  nodeproto.input(0);
   auto it = tensor_type.find(input_name);
   if (it != tensor_type.end()){
      input_type = it->second;
   }else{
      throw std::runtime_error("TMVA::SOFIE ONNX Parser tranpose op has input tensor" + input_name + " but its type is not yet registered");
   }

   std::unique_ptr<ROperator> op;
   std::vector<int_t> attr_perm;

   if (nodeproto.attribute_size() == 1){
      attr_perm.assign(nodeproto.attribute(0).ints().begin(), nodeproto.attribute(0).ints().end());
   }

   switch(input_type){
   case ETensorType::FLOAT:
      if (!attr_perm.empty()){
         op.reset(new ROperator_Transpose<float>(attr_perm, nodeproto.input(0), nodeproto.output(0)));
      }else{
         op.reset(new ROperator_Transpose<float> (nodeproto.input(0), nodeproto.output(0)));
      }
      break;
   default:
      throw std::runtime_error("TMVA::SOFIE - Unsupported - Operator Transpose does not yet support input type " + std::to_string(static_cast<int>(input_type)));
   }

   ETensorType output_type = (op->TypeInference({input_type}))[0];
   auto it2 = tensor_type.find(nodeproto.output(0));
   if (it2 == tensor_type.end()){
      tensor_type[nodeproto.output(0)] = output_type;
   }


   return std::move(op);
}

std::unique_ptr<ROperator> make_ROperator_Relu(const onnx::NodeProto& nodeproto, const onnx::GraphProto& graphproto, std::unordered_map<std::string, ETensorType>& tensor_type){

   ETensorType input_type;

   auto input_name =  nodeproto.input(0);
   auto it = tensor_type.find(input_name);
   if (it != tensor_type.end()){
      input_type = it->second;
   }else{
      throw std::runtime_error("TMVA::SOFIE ONNX Parser relu op has input tensor" + input_name + " but its type is not yet registered");
   }

   std::unique_ptr<ROperator> op;


   switch(input_type){
   case ETensorType::FLOAT:
      op.reset(new ROperator_Relu<float>(nodeproto.input(0), nodeproto.output(0)));
      break;
   default:
      throw std::runtime_error("TMVA::SOFIE - Unsupported - Operator Relu does not yet support input type " + std::to_string(static_cast<int>(input_type)));
   }

   ETensorType output_type = (op->TypeInference({input_type}))[0];
   auto it2 = tensor_type.find(nodeproto.output(0));
   if (it2 == tensor_type.end()){
      tensor_type[nodeproto.output(0)] = output_type;
   }

   return std::move(op);
}

std::unique_ptr<ROperator> make_ROperator_Gemm(const onnx::NodeProto& nodeproto, const onnx::GraphProto& graphproto, std::unordered_map<std::string, ETensorType>& tensor_type){

   ETensorType input_type;

   auto input_name =  nodeproto.input(0);
   auto it = tensor_type.find(input_name);
   if (it != tensor_type.end()){
      input_type = it->second;
   }else{
      throw std::runtime_error("TMVA::SOFIE ONNX Parser gemm op has input tensor" + input_name + " but its type is not yet registered");
   }

   std::unique_ptr<ROperator> op;

   float attr_alpha =1.0;
   float attr_beta =1.0;
   int_t attr_transA =0;
   int_t attr_transB =0;

   for (int i = 0; i < nodeproto.attribute_size(); i++){
      std::string attribute_name = nodeproto.attribute(i).name();
      if (attribute_name == "alpha"){
         attr_alpha = nodeproto.attribute(i).f();
      }else if(attribute_name == "beta"){
         attr_beta = nodeproto.attribute(i).f();
      }else if(attribute_name == "transA"){
         attr_transA = nodeproto.attribute(i).i();
         if (attr_transA != 0 && attr_transA != 1) throw std::runtime_error("TMVA::SOFIE Error - Model Loading - attribute transA in Operator Gemm not 0/1");
      }else if(attribute_name == "transB"){
         attr_transB = nodeproto.attribute(i).i();
         if (attr_transB != 0 && attr_transB != 1) throw std::runtime_error("TMVA::SOFIE Error - Model Loading - attribute transB in Operator Gemm not 0/1");
      }else{
         std::cout << "TMVA::SOFIE Warning - Model Loading - Attribute " << attribute_name << " in OperatorNode " << nodeproto.name() << " is not defined in ONNX IR and not applied!\n";
      }
   }


   switch(input_type){
   case ETensorType::FLOAT:
      if (nodeproto.input_size() == 2){
         op.reset(new ROperator_Gemm<float>(attr_alpha, attr_beta, attr_transA, attr_transB, nodeproto.input(0), nodeproto.input(1), nodeproto.output(0)));
      }else{
         op.reset(new ROperator_Gemm<float>(attr_alpha, attr_beta, attr_transA, attr_transB, nodeproto.input(0), nodeproto.input(1), nodeproto.input(2), nodeproto.output(0)));
      }
      break;
   default:
      throw std::runtime_error("TMVA::SOFIE - Unsupported - Operator Relu does not yet support input type " + std::to_string(static_cast<int>(input_type)));
   }

   ETensorType output_type = (op->TypeInference({input_type, input_type}))[0];
   auto it2 = tensor_type.find(nodeproto.output(0));
   if (it2 == tensor_type.end()){
      tensor_type[nodeproto.output(0)] = output_type;
   }


   return std::move(op);
}

} //INTERNAL



RModel RModelParser_ONNX::Parse(std::string filename){
   char sep = '/';
   #ifdef _WIN32
      sep = '\\';
   #endif
   size_t i = filename.rfind(sep, filename.length());
   std::string modelname;
   if (i != std::string::npos){
      filename = (filename.substr(i+1, filename.length() - i));
   }



   std::time_t ttime = std::time(0);
   std::tm* gmt_time = std::gmtime(&ttime);
   std::string parsetime (std::asctime(gmt_time));




   GOOGLE_PROTOBUF_VERIFY_VERSION;
   //model I/O
   onnx::ModelProto model;
   RModel rmodel(filename, parsetime);

   std::unordered_map<std::string, ETensorType> tensor_type;

   std::fstream input(filename, std::ios::in | std::ios::binary);
   if (!model.ParseFromIstream(&input)){
      throw std::runtime_error("TMVA::SOFIE - Failed to parse onnx file");
   }

   const onnx::GraphProto& graph = model.graph(); //not a memory leak. model freed automatically at the end.
   google::protobuf::ShutdownProtobufLibrary();

   std::unordered_set<std::string> initializer_names;
   for (int i=0; i < graph.initializer_size(); i++){
      initializer_names.insert(graph.initializer(i).name());
   }


   for (int i=0; i < graph.input_size(); i++){

      tensor_type[graph.input(i).name()] = static_cast<ETensorType>(graph.input(i).type().tensor_type().elem_type());

      if (initializer_names.find(graph.input(i).name()) != initializer_names.end())  continue;

      //input datanode is not a weight node (has no initializer)
      const onnx::ValueInfoProto& valueinfoproto = graph.input(i);
      std::string input_name = valueinfoproto.name();
      ETensorType type = static_cast<ETensorType>(valueinfoproto.type().tensor_type().elem_type());
      if (type != ETensorType::FLOAT){
         throw std::runtime_error("TMVA::SOFIE Data type in input tensor " + input_name + " not supported!\n");
      }

      std::vector<Dim> fShape;
      bool existParam = false;
      if (!valueinfoproto.type().tensor_type().has_shape()) throw std::runtime_error("TMVA::SOFIE datanode with no shape restrictions is not supported yet");
      for (int i = 0; i < valueinfoproto.type().tensor_type().shape().dim_size(); i++){
         Dim dim;
         if (valueinfoproto.type().tensor_type().shape().dim(i).value_case() == onnx::TensorShapeProto_Dimension::ValueCase::kDimValue){
            dim.dim = valueinfoproto.type().tensor_type().shape().dim(i).dim_value();
         }else if (valueinfoproto.type().tensor_type().shape().dim(i).value_case() == onnx::TensorShapeProto_Dimension::ValueCase::kDimParam){
            dim.isParam = true;
            existParam = true;
            dim.param = valueinfoproto.type().tensor_type().shape().dim(i).dim_param();
         }else{
            throw std::runtime_error("TMVA::SOFIE ONNX file error: Valueinfoproto " + input_name + " has neither dim_value nor dim_param! \n");
         }
         fShape.push_back(dim);
      }
      if (valueinfoproto.type().tensor_type().shape().dim_size() == 0){
         Dim dim;
         dim.dim = 1;
         fShape.push_back(dim);
      } //in case this TensorShapeProto has no dimension message: ONNX IR defines this to be a scalar

      if (!existParam){
         std::vector<size_t> fShape_sizet;
         for (auto& i: fShape){
            fShape_sizet.push_back(i.dim);
         }

         rmodel.AddInputTensorInfo(input_name, type, fShape_sizet);
      }else{
         rmodel.AddInputTensorInfo(input_name, type, fShape);
      }

   }

   for (int i=0; i < graph.initializer_size(); i++){
      onnx::TensorProto* tensorproto = const_cast<onnx::TensorProto*>(&graph.initializer(i));
      std::vector<std::size_t> fShape;
      std::size_t fLength = 1;
      for (int i = 0; i < tensorproto->dims_size(); i++){
         fShape.push_back(tensorproto->dims(i));
         fLength *= tensorproto->dims(i);
      }

      std::string input_name = graph.initializer(i).name();

      switch(static_cast<ETensorType>(graph.initializer(i).data_type())){
         case ETensorType::FLOAT : {
            //void* data = malloc (fLength * sizeof(float));
            std::shared_ptr<void> data(malloc(fLength * sizeof(float)), free);

            if (tensorproto->raw_data().empty() == false){
               auto raw_data_ptr = reinterpret_cast<float*>(const_cast<char*>(tensorproto->raw_data().c_str()));
               std::memcpy(data.get(), raw_data_ptr, fLength * sizeof(float));
            }else{
               tensorproto->mutable_float_data()->ExtractSubrange(0, tensorproto->float_data_size(), static_cast<float*>(data.get()));
            }

            rmodel.AddInitializedTensor(input_name, ETensorType::FLOAT, fShape, data);
            break;
         }
         default: throw std::runtime_error("Data type in weight tensor " + graph.initializer(i).name() + " not supported!\n");
      }
   }



   for (int i=0; i < graph.node_size(); i++){
      rmodel.AddOperator(std::move(INTERNAL::make_ROperator(i, graph, tensor_type)));
   }

   std::vector<std::string> outputnames;
   for (int i=0; i < graph.output_size(); i++){
      outputnames.push_back(graph.output(i).name());
   }
   rmodel.AddOutputTensorNameList(outputnames);

   return rmodel;

}



}//SOFIE
}//Experimental
}//TMVA
