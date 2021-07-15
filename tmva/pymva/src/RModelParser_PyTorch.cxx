#include "TMVA/RModelParser_PyTorch.h"

namespace TMVA{
namespace Experimental{
namespace SOFIE{

std::unordered_map<std::string, NodeType> NType =
    {
        {"'onnx::gemm'", NodeType::GEMM},
        {"'onnx::relu'", NodeType::RELU},
        {"'onnx::transpose'", NodeType::TRANSPOSE}
    };




namespace PyTorch{

RModel Parse(std::string filename, std::vector<std::vector<size_t>> inputShapes, std::vector<ETensorType> inputDTypes){
    char sep='/';
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


    //Extracting model information
    //Model is converted to ONNX graph format
    //using PyTorch's internal function with the input shape provided
    PyRunString("globals().update(locals())",fGlobalNS,fLocalNS);
    PyRunString("import torch",fGlobalNS,fLocalNS);
    PyRunString("print('Torch Version: '+torch.__version__)",fGlobalNS,fLocalNS);
    PyRunString("from torch.onnx.utils import _model_to_graph",fGlobalNS,fLocalNS);
    PyRunString(TString::Format("model= torch.jit.load('%s')",filename.c_str()),fGlobalNS,fLocalNS);
    PyRunString("globals().update(locals())",fGlobalNS,fLocalNS);
    PyRunString("model.eval()",fGlobalNS,fLocalNS);

    //Building dummy inputs for the model
    PyRunString("dummy=[]",fGlobalNS,fLocalNS);
    for(auto it=0;it<inputShapes.size();++it){
        PyRunString("inputShape=[]",fGlobalNS,fLocalNS);
        for(auto itr=0;itr<inputShapes[it].size();++itr){
            PyRunString(TString::Format("inputShape.append(%d)",(int)inputShapes[it][itr]),fGlobalNS,fLocalNS);
            }
        switch(inputDTypes[it]){
            case ETensorType::FLOAT:{
                PyRunString("dummy.append(torch.rand(*inputShape))",fGlobalNS,fLocalNS);
                break;
            }
            default:
                throw std::runtime_error("Type Error: TMVA SOFIE does not yet support the input tensor data type"+ConvertTypeToString(inputDTypes[it]));
        }
        }

    //Finding example outputs from dummy
    PyRunString("output=model(*dummy)",fGlobalNS,fLocalNS);

    //Getting the ONNX graph from model using the dummy inputs and example outputs
    PyRunString("graph=_model_to_graph(model,dummy,example_outputs=output)",fGlobalNS,fLocalNS);



    //Extracting Input tensor info
    PyRunString("inputs=[x for x in model.graph.inputs()]",fGlobalNS,fLocalNS);
    PyRunString("inputs=inputs[1:]",fGlobalNS,fLocalNS);
    PyRunString("inputNames=[x.debugName() for x in inputs]",fGlobalNS,fLocalNS);
    PyObject* pInputs= PyDict_GetItemString(fLocalNS,"inputNames");
    for(Py_ssize_t inputIter=0; inputIter<PyList_Size(pInputs);++inputIter){
        std::string inputName(PyStringAsString(PyList_GetItem(pInputs,inputIter)));
        std::vector<size_t>inputShape=inputShapes[inputIter];
        ETensorType inputDType=inputDTypes[inputIter];
        switch(inputDType){
            case(ETensorType::FLOAT): {
                rmodel.AddInputTensorInfo(inputName, ETensorType::FLOAT, inputShape);
                break;
            }
            default:
                throw std::runtime_error("Type Error: TMVA SOFIE does not yet support the input tensor data type"+ConvertTypeToString(inputDType));
        }
        }
    Py_XDECREF(pInputs);



    //Extracting the model information in list modelData
    PyRunString("modelData=[]",fGlobalNS,fLocalNS);
    PyRunString("for i in graph[0].nodes():\n"
                "    globals().update(locals())\n"
                "    nodeData=[]\n"
                "    nodeData.append(i.kind())\n"
                "    nodeAttributeNames=[x for x in i.attributeNames()]\n"
                "    nodeAttributes={j:i[j] for j in nodeAttributeNames}\n"
                "    nodeData.append(nodeAttributes)\n"
                "    nodeInputs=[x for x in i.inputs()]\n"
                "    nodeInputNames=[x.debugName() for x in nodeInputs]\n"
                "    nodeData.append(nodeInputNames)\n"
                "    nodeOutputs=[x for x in i.outputs()]\n"
                "    nodeOutputNames=[x.debugName() for x in nodeOutputs]\n"
                "    nodeData.append(nodeOutputNames)\n"
                "    nodeDType=[x.type().scalarType() for x in nodeOutputs]\n"
                "    nodeData.append(nodeDType)\n"
                "    modelData.append(nodeData)",fGlobalNS,fLocalNS);
    Py_ssize_t modelIterator, modelSize;
    PyObject* pModel =PyDict_GetItemString(fLocalNS,"modelData");
    PyObject* node,*attributes,*inputs,*outputs,*nodeDType;
    modelSize = PyList_Size(pModel);

    for(modelIterator=0;modelIterator<modelSize;++modelIterator){
        node=PyList_GetItem(pModel,modelIterator);
        std::string type(PyStringAsString(PyList_GetItem(node,0)));

        if(NType.find(toLower(type))==NType.end())
            throw std::runtime_error("Layer error: TMVA SOFIE does not yet suppport layer type"+type);

        attributes=PyList_GetItem(node,1);
        inputs=PyList_GetItem(node,2);
        outputs=PyList_GetItem(node,3);
        ETensorType nodeDType=convertStringToType(PyStringAsString(PyList_GetItem(PyList_GetItem(node,4),0)));

        switch(NType.find(toLower(type))->second){
            case NodeType::GEMM : {
                float attr_alpha = (float)(PyFloat_AsDouble(PyDict_GetItemString(attributes,"alpha")));
                float attr_beta = (float)(PyFloat_AsDouble(PyDict_GetItemString(attributes,"beta")));
                int_t attr_transA;
                int_t attr_transB;

                if(PyDict_Contains(attributes,PyUnicode_FromString("transB"))){
                         attr_transB=(int_t)(PyLong_AsLong(PyDict_GetItemString(attributes,"transB")));
                         attr_transA=!attr_transB;
                    }
                else{
                        attr_transA=(int_t)(PyLong_AsLong(PyDict_GetItemString(attributes,"transA")));
                        attr_transB=!attr_transA;
                    }
                switch(nodeDType){
                case ETensorType::FLOAT: {
                    std::unique_ptr<ROperator> op;
                    op.reset(new ROperator_Gemm<float>(attr_alpha, attr_beta, attr_transA, attr_transB, PyStringAsString(PyList_GetItem(inputs,0)), PyStringAsString(PyList_GetItem(inputs,1)), PyStringAsString(PyList_GetItem(inputs,2)), PyStringAsString(PyList_GetItem(outputs,0))));
                    rmodel.AddOperator(std::move(op));
                    break;
                    }
                default:
                    throw std::runtime_error("TMVA::SOFIE - Unsupported - Operator Gemm does not yet support input type " + ConvertTypeToString(nodeDType));

                }
                break;
                }

            case NodeType::RELU : {
                switch(nodeDType){
                    case ETensorType::FLOAT: {
                        std::unique_ptr<ROperator> op;
                        op.reset(new ROperator_Relu<float>(PyStringAsString(PyList_GetItem(inputs,0)), PyStringAsString(PyList_GetItem(outputs,0))));
                        rmodel.AddOperator(std::move(op));
                        break;
                    }
                    default:
                    throw std::runtime_error("TMVA::SOFIE - Unsupported - Operator Relu does not yet support input type " + ConvertTypeToString(nodeDType));
                }
                break;
                }

            case NodeType::TRANSPOSE:{
                std::vector<int_t> attr_perm;
                PyObject* permute=PyDict_GetItemString(attributes,"perm");
                for(Py_ssize_t permIter=0; permIter<PyList_Size(permute);++permIter){
                    attr_perm.push_back((int_t)PyLong_AsLong(PyList_GetItem(permute,permIter)));
                }
                switch(nodeDType){
                case ETensorType::FLOAT: {
                std::unique_ptr<ROperator> op;
                op.reset(new ROperator_Transpose<float>(attr_perm, PyStringAsString(PyList_GetItem(inputs,0)), PyStringAsString(PyList_GetItem(outputs,0))));
                rmodel.AddOperator(std::move(op));
                break;
                }
                default:
                    throw std::runtime_error("TMVA::SOFIE - Unsupported - Operator Transpose does not yet support input type " + ConvertTypeToString(nodeDType));
                }

                Py_XDECREF(permute);
                break;
                }

            default:
                throw std::runtime_error("Node Error: TMVA SOFIE does not yet support node type " + type);
            }
            }
    Py_XDECREF(nodeDType);
    Py_XDECREF(outputs);
    Py_XDECREF(inputs);
    Py_XDECREF(attributes);
    Py_XDECREF(node);
    Py_XDECREF(pModel);


    //Extracting model weights to add the initialized tensors to the RModel
    PyRunString("weightNames=[k for k in graph[1].keys()]",fGlobalNS,fLocalNS);
    PyRunString("weights=[v.numpy() for v in graph[1].values()]",fGlobalNS,fLocalNS);
    PyRunString("weightDTypes=[v.numpy().dtype.name for v in graph[1].values()]",fGlobalNS,fLocalNS);
    PyObject* weightNames= PyDict_GetItemString(fLocalNS,"weightNames");
    PyObject* weightTensors= PyDict_GetItemString(fLocalNS,"weights");
    PyObject* weightDTypes= PyDict_GetItemString(fLocalNS,"weightDTypes");
    PyObject* weightTensor;
    for(Py_ssize_t weightIter=0; weightIter<PyList_Size(weightNames);++weightIter){
        weightTensor= PyList_GetItem(weightTensors,weightIter);
        std::string weightName(PyStringAsString(PyList_GetItem(weightNames,weightIter)));
        ETensorType weightDType=convertStringToType(PyStringAsString(PyList_GetItem(weightDTypes,weightIter)));

        switch(weightDType){
            case ETensorType::FLOAT:{
                //Converting the numpy array object to RTensor
                RTensor<float> value=getArray(weightTensor);
                std::shared_ptr<void> data(malloc(value.GetSize() * sizeof(float)), free);
                std::memcpy(data.get(),value.GetData(),value.GetSize() * sizeof(float));
                rmodel.AddInitializedTensor(weightName, ETensorType::FLOAT,value.GetShape(), data);
                break;
                }
            default:
                throw std::runtime_error("Type error: TMVA SOFIE does not yet weight data layer type"+ConvertTypeToString(weightDType));
            }
    }

    Py_XDECREF(weightDTypes);
    Py_XDECREF(weightNames);
    Py_XDECREF(weightTensors);
    Py_XDECREF(weightTensor);


    //Extracting output tensor names
    PyRunString("outputs=[x for x in graph[0].outputs()]",fGlobalNS,fLocalNS);
    PyRunString("outputNames=[x.debugName() for x in outputs]",fGlobalNS,fLocalNS);
    PyObject* pOutputs= PyDict_GetItemString(fLocalNS,"outputNames");
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

RModel Parse(std::string filepath,std::vector<std::vector<size_t>> inputShapes){
      std::vector<ETensorType> dtype(inputShapes.size(),ETensorType::FLOAT);
      return Parse(filepath,inputShapes,dtype);
    }
}
}
}
}
