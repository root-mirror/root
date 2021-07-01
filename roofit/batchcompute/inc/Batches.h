// Author: Emmanouil Michalainas, CERN 3 March 2021

#ifndef ROOFIT_BATCHCOMPUTE_BATCHES_H
#define ROOFIT_BATCHCOMPUTE_BATCHES_H

#include "rbc.h"

#include <stdint.h>

namespace RooBatchCompute {

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
    
    __device__ constexpr bool isItVector() { return isVector; }
    inline void set(double _scalar, InputArr _array, bool _isVector)
    {
      scalar = _scalar;
      array = _array;
      isVector = _isVector;
    }
    inline void advance(size_t _nEvents) { array += isVector*_nEvents; }
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
    __device__ constexpr size_t getNEvents() const { return nEvents; }
    __device__ constexpr uint8_t getNExtraArgs() const { return nExtraArgs; }
    __device__ constexpr double extraArg(uint8_t i) const { return extraArgs[i]; }
    __device__ constexpr Batch operator[] (int batchIdx) const { return arrays[batchIdx]; }
    inline void setNEvents(size_t n=bufferSize) { nEvents = n; }
    inline void advance(size_t _nEvents)
    {
      for (int i=0; i<nBatches; i++) arrays[i].advance(_nEvents);
      output += _nEvents;
    }
}; //end class Batches
} // End namespace RF_ARCH
} //end namespace RooBatchCompute
#endif // #ifdef ROOFIT_BATCHCOMPUTE_BATCHES_H 
