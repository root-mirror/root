// RooBatchCompute library created September 2020 by Emmanouil Michalainas

#include "RooBatchCompute.h"
#include "Batches.h"

#include "TEnv.h"
#include "TError.h"

#include <thrust/reduce.h>

namespace RooBatchCompute {
namespace RF_ARCH {
  
std::vector<void(*)(Batches)> getFunctions();

class RooBatchComputeClass : public RooBatchComputeInterface {
  private:
    const std::vector<void(*)(Batches)> computeFunctions;
  public:
    RooBatchComputeClass()
      : computeFunctions(getFunctions())
    {
      dispatch_gpu = this; // Set the dispatch pointer to this instance of the library upon loading
    }

    void init()
    {
      cudaError_t err = cudaSetDevice(0);
      if (err==cudaSuccess) cudaFree(nullptr);
      else
      {
        dispatch_gpu = nullptr;
        Error( (std::string(__func__)+"(), "+__FILE__+":"+std::to_string(__LINE__)).c_str(), "%s", cudaGetErrorString(err) );
      }
    }

    void compute(Computer computer, RestrictArr output, size_t nEvents, const DataMap& varData, const VarVector& vars, const ArgVector& extraArgs) override  
    {
      Batches batches(output, nEvents, varData, vars, extraArgs);
      computeFunctions[computer]<<<256,64>>>(batches);
    }

    double sumReduce(InputArr input, size_t n) override
    {
      return thrust::reduce(thrust::device, input, input+n, 0.0);
    }

    void* malloc(size_t size) override
    {
      void* ret = nullptr;
      cudaError_t error = cudaMalloc(&ret, size);
      if (error != cudaSuccess)
      {
        Fatal( (std::string(__func__)+"(), "+__FILE__+":"+std::to_string(__LINE__)).c_str(), "%s", cudaGetErrorString(error) );
        throw std::bad_alloc();
      }
      else return ret;
    }

    void free(void* ptr) override
    {
      cudaFree(ptr);
    }

    void memcpyToGPU(void* dest, const void* src, size_t n)
    {
      cudaMemcpy(dest, src, n, cudaMemcpyHostToDevice);
    }
    
    void memcpyToCPU(void* dest, const void* src, size_t n)
    {
      cudaMemcpy(dest, src, n, cudaMemcpyDeviceToHost);
    }
}; // End class RooBatchComputeClass

/// Static object to trigger the constructor which overwrites the dispatch pointer.
static RooBatchComputeClass computeObj;


Batches::Batches(RestrictArr _output, size_t _nEvents, const DataMap& varData, const VarVector& vars, const ArgVector& _extraArgs, double[maxParams][bufferSize])
  : nEvents(_nEvents), nBatches(vars.size()), nExtraArgs(_extraArgs.size()), output(_output)
{  
  for (int i=0; i<vars.size(); i++)
  {
    const RooSpan<const double>& span = varData.at(vars[i]);
    size_t size = span.size();
    if (size==1) arrays[i].set(span[0], nullptr, false);
    else arrays[i].set(0.0, span.data(), true);
  }
  std::copy(_extraArgs.cbegin(), _extraArgs.cend(), extraArgs);
}

} // End namespace RF_ARCH
} // End namespace RooBatchCompute
