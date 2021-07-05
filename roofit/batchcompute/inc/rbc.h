#ifndef RBC_H
#define RBC_H

#include "RooSpan.h"
#include <unordered_map>
#include <vector>

#ifndef __CUDACC__
  #define __device__
  #define __global__
  #define __host__
#endif // #ifndef __CUDACC__

class RooAbsReal;

namespace RooBatchCompute {

struct RunContext;
typedef std::unordered_map<const RooAbsReal*,RooSpan<const double>> DataMap;
typedef std::vector<const RooAbsReal*> VarVector;
typedef std::vector<double> ArgVector;
typedef double* __restrict RestrictArr;
typedef const double* __restrict InputArr;

}
namespace rbc=RooBatchCompute;

#endif //#ifndef RBC_H
