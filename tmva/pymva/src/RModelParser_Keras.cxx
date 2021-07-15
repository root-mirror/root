#include "TMVA/RModelParser_Keras.h"

namespace TMVA{
namespace Experimental{
namespace SOFIE{


std::unordered_map<std::string, LayerType> Type =
    {
        {"'dense'", LayerType::DENSE},
        {"'relu'", LayerType::RELU},
        {"'permute'", LayerType::TRANSPOSE}
    };


namespace INTERNAL{

   std::unique_ptr<ROperator> make_ROperator_Gemm(std::string input,std::string output,std::string kernel,std::string bias, ETensorType dtype)
   {
      std::unique_ptr<ROperator> op;

      float attr_alpha =1.0;
      float attr_beta =1.0;
      int_t attr_transA =1;
      int_t attr_transB =0;

      switch(dtype){
         case ETensorType::FLOAT:
         op.reset(new ROperator_Gemm<float>(attr_alpha, attr_beta, attr_transA, attr_transB, input, kernel, bias, output));
         break;

         default:
         throw std::runtime_error("TMVA::SOFIE - Unsupported - Operator Gemm does not yet support input type " + ConvertTypeToString(dtype));
         }

         return std::move(op);
   }

   std::unique_ptr<ROperator> make_ROperator_Relu(std::string input, std::string output, ETensorType dtype)
   {
      std::unique_ptr<ROperator> op;
      switch(dtype){
         case ETensorType::FLOAT:
         op.reset(new ROperator_Relu<float>(input, output));
         break;
         default:
         throw std::runtime_error("TMVA::SOFIE - Unsupported - Operator Relu does not yet support input type " + ConvertTypeToString(dtype));
         }
   return std::move(op);
   }

   std::unique_ptr<ROperator> make_ROperator_Transpose(std::string input, std::string output, std::vector<int_t> dims, ETensorType dtype)
   {
      std::unique_ptr<ROperator> op;
      std::vector<int_t> attr_perm=dims;
      switch(dtype){
         case ETensorType::FLOAT:
            if (!attr_perm.empty()){
               op.reset(new ROperator_Transpose<float>(attr_perm, input, output));
               }
            else{
               op.reset(new ROperator_Transpose<float> (input, output));
               }
         break;
         default:
            throw std::runtime_error("TMVA::SOFIE - Unsupported - Operator Transpose does not yet support input type " + ConvertTypeToString(dtype));
            }
   return std::move(op);
   }

}



namespace PyKeras {

RModel Parse(std::string filename){

   char sep = '/';
   #ifdef _WIN32
   sep = '\\';
   #endif

   size_t i = filename.rfind(sep, filename.length());
   if (i != std::string::npos){
      filename = (filename.substr(i+1, filename.length() - i));
   }

   if(!std::ifstream(filename).good()){
        throw std::runtime_error("Model file "+filename+" not found!");
    }


   std::time_t ttime = std::time(0);
   std::tm* gmt_time = std::gmtime(&ttime);
   std::string parsetime (std::asctime(gmt_time));

   RModel rmodel(filename, parsetime);


   Py_Initialize();

   PyObject* main = PyImport_AddModule("__main__");
   PyObject* fGlobalNS = PyModule_GetDict(main);
   PyObject* fLocalNS = PyDict_New();
   if (!fGlobalNS) {
       throw std::runtime_error("Can't init global namespace for Python");
       }
   if (!fLocalNS) {
       throw std::runtime_error("Can't init local namespace for Python");
       }

   //Extracting model information: For each layer: type,name,activation,dtype,input tensor's name,
   //output tensor's name, kernel's name, bias's name
   //None object is returned for if property doesn't belong to layer
   PyRunString("import keras",fGlobalNS,fLocalNS);
   PyRunString("from keras.models import load_model",fGlobalNS,fLocalNS);
   PyRunString("print('Keras Version: '+ keras.__version__)",fGlobalNS,fLocalNS);
   PyRunString(TString::Format("model=load_model('%s')",filename.c_str()),fGlobalNS,fLocalNS);
   PyRunString(TString::Format("model.load_weights('%s')",filename.c_str()),fGlobalNS,fLocalNS);
   PyRunString("globals().update(locals())",fGlobalNS,fLocalNS);
   PyRunString("modelData=[]",fGlobalNS,fLocalNS);
   PyRunString("for idx in range(len(model.layers)):\n"
            "  layer=model.get_layer(index=idx)\n"
            "  globals().update(locals())\n"
            "  layerData=[]\n"
            "  layerData.append(layer.__class__.__name__)\n"
            "  layerData.append(layer.get_config())\n"
            "  layerData.append(layer.input)\n"
            "  layerData.append(layer.output)\n"
            "  layerData.append([x.name for x in layer.weights])\n"
            "  modelData.append(layerData)",fGlobalNS,fLocalNS);


   Py_ssize_t modelIterator, modelSize;
   PyObject* pModel = PyDict_GetItemString(fLocalNS,"modelData");
   PyObject* layer,*attributes,*inputs,*outputs;
   modelSize = PyList_Size(pModel);

   for(modelIterator=0;modelIterator<modelSize;++modelIterator){
      layer=PyList_GetItem(pModel,modelIterator);

      std::string type(PyStringAsString(PyList_GetItem(layer,0)));
      attributes=PyList_GetItem(layer,1);
      inputs=PyList_GetItem(layer,2);
      outputs=PyList_GetItem(layer,3);
      ETensorType dtype = convertStringToType(PyStringAsString(PyDict_GetItemString(attributes,"dtype")));

      switch(Type.find(toLower(type))->second){
         case LayerType::DENSE : {

         std::string activation(PyStringAsString(PyDict_GetItemString(attributes,"activation")));
         std::string name(PyStringAsString(PyDict_GetItemString(attributes,"name")));
         std::string input(PyStringAsString(PyObject_GetAttrString(inputs,"name")));
         std::string output(PyStringAsString(PyObject_GetAttrString(outputs,"name")));

         PyObject* weightNames=PyList_GetItem(layer,4);
         std::string kernel(PyStringAsString(PyList_GetItem(weightNames,0)));
         std::string bias(PyStringAsString(PyList_GetItem(weightNames,1)));

                  if(activation != "'linear'"){
                     rmodel.AddOperator(std::move(INTERNAL::make_ROperator_Gemm(input,name+"_gemm",kernel,bias,dtype)));

                     if(Type.find(toLower(activation))==Type.end())
                       throw std::runtime_error("Type error: Layer activation type "+toLower(activation)+" not yet registered in TMVA SOFIE");

                     switch(Type.find(toLower(activation))->second){
                        case LayerType::RELU: {
                           rmodel.AddOperator(std::move(INTERNAL::make_ROperator_Relu(name+"_gemm",output,dtype)));
                           break;
                        }
                        default: throw std::runtime_error("Activation error: TMVA SOFIE does not yet suppport Activation type"+activation);
                        }
                        }
                  else{
                     rmodel.AddOperator(std::move(INTERNAL::make_ROperator_Gemm(input,output,kernel,bias,dtype)));
                  }

                  Py_XDECREF(weightNames);
                  break;
               }

         case LayerType::RELU: {
            std::string input(PyStringAsString(PyObject_GetAttrString(inputs,"name")));
            std::string output(PyStringAsString(PyObject_GetAttrString(outputs,"name")));
            rmodel.AddOperator(std::move(INTERNAL::make_ROperator_Relu(input,output,dtype)));  break;
            }
         case LayerType::TRANSPOSE: {
            std::string input(PyStringAsString(PyObject_GetAttrString(inputs,"name")));
            std::string output(PyStringAsString(PyObject_GetAttrString(outputs,"name")));
            PyObject* permute=PyDict_GetItemString(attributes,"dims");
            std::vector<int_t>dims;
            for(Py_ssize_t tupleIter=0;tupleIter<PyTuple_Size(permute);++tupleIter)
               dims.push_back((int_t)PyLong_AsLong(PyTuple_GetItem(permute,tupleIter)));
            rmodel.AddOperator(std::move(INTERNAL::make_ROperator_Transpose(input,output,dims,dtype)));

            Py_XDECREF(permute);
            break;
            }
         default: throw std::runtime_error("Layer error: TMVA SOFIE does not yet suppport layer type"+type);
         }
         }

   Py_XDECREF(outputs);
   Py_XDECREF(inputs);
   Py_XDECREF(attributes);
   Py_XDECREF(layer);
   Py_XDECREF(pModel);


   //Extracting model's weights
   //For every initialized tensor, weightProp will have its name and dtype in string
   //and value in numpy array
   PyRunString("weight=[]",fGlobalNS,fLocalNS);
   PyRunString("for idx in range(len(model.get_weights())):\n"
               "  weightProp={}\n"
               "  weightProp['name']=model.weights[idx].name\n"
               "  weightProp['dtype']=(model.get_weights())[idx].dtype.name\n"
               "  weightProp['value']=(model.get_weights())[idx]\n"
               "  weight.append(weightProp)",fGlobalNS,fLocalNS);

   PyObject *weightTensor,*weightValue;
   PyObject* pWeight = PyDict_GetItemString(fLocalNS,"weight");

   for (Py_ssize_t weightIter = 0; weightIter < PyList_Size(pWeight); weightIter++) {
      weightTensor  = PyList_GetItem(pWeight, weightIter);
      std::string weightName(PyStringAsString(PyDict_GetItemString(weightTensor,"name")));
      ETensorType weightType= convertStringToType(PyStringAsString(PyDict_GetItemString(weightTensor,"dtype")));
      weightValue   = PyDict_GetItemString(weightTensor,"value");

      //Converting numpy array to RTensor
      RTensor<float> value = getArray(weightValue);

   switch(weightType){
       case ETensorType::FLOAT : {
       std::shared_ptr<void> data(malloc(value.GetSize() * sizeof(float)), free);
       std::memcpy(data.get(),value.GetData(),value.GetSize() * sizeof(float));
       rmodel.AddInitializedTensor(weightName, ETensorType::FLOAT,value.GetShape(), data);
       break;
       }
       default:
          throw std::runtime_error("Type error: TMVA SOFIE does not yet weight data layer type"+ConvertTypeToString(weightType));
      }
     }

   Py_XDECREF(weightTensor);
   Py_XDECREF(weightValue);
   Py_XDECREF(pWeight);


   //Extracting input tensor info
   //For every input tensor inputNames will have their names as string,inputShapes will have their
   //shape as Python Tuple, and inputTypes will have their dtype as string
   PyRunString("inputNames=model.input_names",fGlobalNS,fLocalNS);
   PyRunString("inputShapes=model.input_shape",fGlobalNS,fLocalNS);
   PyRunString("inputTypes=[]",fGlobalNS,fLocalNS);
   PyRunString("for idx in range(len(model.inputs)):\n"
               "  inputTypes.append(model.inputs[idx].dtype.__str__()[9:-2])",fGlobalNS,fLocalNS);

   PyObject* pInputs   = PyDict_GetItemString(fLocalNS,"inputNames");
   PyObject* pInputShapes  = PyDict_GetItemString(fLocalNS,"inputShapes");
   PyObject* pInputTypes   = PyDict_GetItemString(fLocalNS,"inputTypes");

   //For single input models, the model.input_shape will return a tuple describing the input tensor shape
   //For multiple inputs models, the model.input_shape will return a list of tuple, each describing the input tensor shape.
   if(PyTuple_Check(pInputShapes)){
      std::string inputName(PyStringAsString(PyList_GetItem(pInputs,0)));
      ETensorType inputDType = convertStringToType(PyStringAsString(PyList_GetItem(pInputTypes,0)));


      switch(inputDType){

         case ETensorType::FLOAT : {
         if (PyTuple_GetItem(pInputShapes,0) == Py_None){
            throw std::runtime_error("None error: Models not initialized with batch-size are not yet supported in TMVA SOFIE");
         }

         std::vector<size_t>inputShape=getShape(pInputShapes);
         rmodel.AddInputTensorInfo(inputName, ETensorType::FLOAT, inputShape);
         break;
         }

         default:
         throw std::runtime_error("Type error: TMVA SOFIE does not yet suppport data type"+ConvertTypeToString(inputDType));

      }

   }

   else{

      for(Py_ssize_t inputIter = 0; inputIter < PyList_Size(pInputs);++inputIter){

      std::string inputName(PyStringAsString(PyList_GetItem(pInputs,inputIter)));
      ETensorType inputDType = convertStringToType(PyStringAsString(PyList_GetItem(pInputTypes,inputIter)));

      switch(inputDType){

         case ETensorType::FLOAT : {
         PyObject* shapeTuple=PyList_GetItem(pInputShapes,inputIter);

         if (PyTuple_GetItem(shapeTuple,0) == Py_None){
            throw std::runtime_error("None error: Models not initialized with batch-size are not yet supported in TMVA SOFIE");
         }

         std::vector<size_t>inputShape=getShape(shapeTuple);
         rmodel.AddInputTensorInfo(inputName, ETensorType::FLOAT, inputShape);
         break;
         }

         default:
         throw std::runtime_error("Type error: TMVA SOFIE does not yet suppport data type"+ConvertTypeToString(inputDType));

      }
      }

   }

   Py_XDECREF(pInputs);
   Py_XDECREF(pInputShapes);
   Py_XDECREF(pInputTypes);

   //Extracting Output Tensor Names
   PyRunString("outputNames=[]",fGlobalNS,fLocalNS);
   PyRunString("for layerName in model.output_names:\n"
               "    outputNames.append(model.get_layer(layerName).output.name)",fGlobalNS,fLocalNS);
   PyObject* pOutputs   = PyDict_GetItemString(fLocalNS,"outputNames");
   std::vector<std::string> outputNames;
   for(Py_ssize_t outputIter = 0; outputIter < PyList_Size(pOutputs);++outputIter){
         outputNames.push_back(PyStringAsString(PyList_GetItem(pOutputs,outputIter)));
   }
   rmodel.AddOutputTensorNameList(outputNames);

   Py_XDECREF(pOutputs);
   Py_XDECREF(fLocalNS);
   Py_XDECREF(fGlobalNS);
   Py_XDECREF(main);
   Py_Finalize();
   return rmodel;

     }
   }
}
}
}
