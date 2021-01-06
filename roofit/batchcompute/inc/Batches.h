// Author: Emmanouil Michalainas, CERN 3 March 2021

#ifndef ROOFIT_BATCHCOMPUTE_BATCHES_H
#define ROOFIT_BATCHCOMPUTE_BATCHES_H

#include "RooSpan.h"

#include <stdint.h>
#include <vector>
#include <unordered_map>

#ifndef __CUDACC__
  #define __device__
  #define __global__
  #define __host__
#endif // #ifndef __CUDACC__
class RooAbsReal;

namespace RooBatchCompute {

typedef std::unordered_map<const RooAbsReal*,RooSpan<const double>> DataMap;
typedef std::vector<const RooAbsReal*> VarVector;
typedef std::vector<double> ArgVector;
typedef double* __restrict RestrictArr;
typedef const double* __restrict InputArr;

constexpr uint8_t maxParams=8;
constexpr uint8_t maxExtraArgs=16;
constexpr uint16_t bufferSize=64;

namespace RF_ARCH {

class Batch {
  private:
    double scalar=0;
    const double* __restrict array=nullptr;
    bool isVector=false;

  public:
    Batch() = default;
    inline Batch(InputArr _array, bool _isVector)
      : scalar{_array[0]}, array{_array}, isVector{_isVector}
    {}

    inline void set(double _scalar, InputArr _array, bool _isVector)
    {
      scalar = _scalar;
      array = _array;
      isVector = _isVector;
    }
    inline void advance() { array += isVector*bufferSize; }
#ifdef __CUDACC__
    __device__ constexpr double operator[](size_t i) noexcept { return isVector ? array[i] : scalar; }
#else
    constexpr double operator[](size_t i) noexcept { return array[i]; }
#endif // #ifdef __CUDACC__
}; //end class Batch

/////////////////////////////////////////////////////////////////////////////////////////////////////////

class Batches {
  private:
    Batch arrays[maxParams];
    size_t nEvents=0;
    double extraArgs[maxExtraArgs];
    uint8_t nBatches=0;
    uint8_t nExtraArgs;

  public:
    RestrictArr output=nullptr;
    
    Batches(RestrictArr output, size_t nEvents, const DataMap& varData, const VarVector& vars, const ArgVector& extraArgs={}, double stackArr[maxParams][bufferSize]=nullptr);
    __device__ constexpr size_t getNEvents() { return nEvents; }
    __device__ constexpr uint8_t getNExtraArgs() { return nExtraArgs; }
    __device__ constexpr double extraArg(uint8_t i) { return extraArgs[i]; }
    __device__ constexpr Batch operator[] (int batchIdx) { return arrays[batchIdx]; }
    inline void setNEvents(size_t n=bufferSize) { nEvents = n; }
    inline void advance()
    {
      for (int i=0; i<nBatches; i++) arrays[i].advance();
      output += bufferSize;
    }
}; //end class Batches
} // End namespace RF_ARCH
} //end namespace RooBatchCompute
#endif // #ifdef ROOFIT_BATCHCOMPUTE_BATCHES_H 
