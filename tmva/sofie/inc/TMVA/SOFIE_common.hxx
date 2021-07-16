#ifndef TMVA_SOFIE_SOFIE_COMMON
#define TMVA_SOFIE_SOFIE_COMMON

#include "TMVA/RTensor.hxx"

#include <type_traits>
#include <cstdint>
#include <string>
#include <vector>
#include <memory>
#include <regex>

namespace TMVA{
namespace Experimental{
namespace SOFIE{

//typedef RTensor tensor_t;

enum class ETensorType{
   UNDEFINED = 0, FLOAT = 1, UNINT8 = 2, INT8 = 3, UINT16 = 4, INT16 = 5, INT32 = 6, INT64 = 7, STRING = 8, BOOL = 9, //order sensitive
    FLOAT16 = 10, DOUBLE = 11, UINT32 = 12, UINT64 = 13, COMPLEX64 = 14, COMPLEX28 = 15, BFLOAT16 = 16
};


typedef std::int64_t int_t;

std::string ConvertTypeToString(ETensorType type);

struct Dim{
   bool isParam = false;
   size_t dim;
   std::string param;
};

std::vector<Dim> ConvertShapeToDim(std::vector<size_t> shape);


struct InputTensorInfo{
   ETensorType type;
   std::vector<Dim> shape;
};

struct TensorInfo{
   ETensorType type;
   std::vector<size_t> shape;
};

std::size_t ConvertShapeToLength(std::vector<size_t> shape);

struct InitializedTensor{
   ETensorType type;
   std::vector<std::size_t> shape;
   std::shared_ptr<void> data;
   //void* data;
};

template <typename T>
ETensorType GetTemplatedType(T obj){
   if (std::is_same<T, float>::value) return ETensorType::FLOAT;
   if (std::is_same<T, uint8_t>::value) return ETensorType::UNINT8;
   if (std::is_same<T, int8_t>::value) return ETensorType::INT8;
   if (std::is_same<T, uint16_t>::value) return ETensorType::UINT16;
   if (std::is_same<T, int16_t>::value) return ETensorType::INT16;
   if (std::is_same<T, int32_t>::value) return ETensorType::INT32;
   if (std::is_same<T, int64_t>::value) return ETensorType::INT64;
   if (std::is_same<T, std::string>::value) return ETensorType::STRING;
   if (std::is_same<T, bool>::value) return ETensorType::BOOL;
   //float16 unimplemented
   if (std::is_same<T, double>::value) return ETensorType::DOUBLE;
   if (std::is_same<T, uint32_t>::value) return ETensorType::UINT32;
   if (std::is_same<T, uint64_t>::value) return ETensorType::UINT64;
   //complex 64, 28, bfloat 16 unimplemented
}

namespace UTILITY{
template<typename T>
T* Unidirectional_broadcast(const T* original_data, const std::vector<size_t> original_shape, const std::vector<size_t> target_shape);
std::string Clean_name(std::string input_tensor_name);
}

namespace BLAS{
extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
                       const float * beta, float * C, const int * ldc);
}//BLAS
}//SOFIE
}//Experimental
}//TMVA

#endif //TMVA_SOFIE_RMODEL
