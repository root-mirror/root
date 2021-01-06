#ifndef ROO_FIT_DRIVER_H
#define ROO_FIT_DRIVER_H

#include "RooBatchCompute.h"

#include <queue>
#include <unordered_map>

class RooAbsData;
class RooAbsReal;
class RooNLLVarNew;

class RooFitDriver {
  public:
     RooFitDriver(const RooAbsData& data, const RooNLLVarNew& topNode, int batchMode);
     ~RooFitDriver();
     //~  inline RooAbsReal* getTopNode() { return initialQueue.back(); }
     double getVal();
     RooArgSet* getParameters();
     
  private:
    const int batchMode=0;
    double* cudaMemDataset=nullptr;
    // used for preserving static info about the computation graph
    rbc::DataMap dataMap;
    const RooArgSet* observables=nullptr;
    const RooNLLVarNew& topNode;
    size_t nEvents;
    std::queue<const RooAbsReal*> initialQueue;
    std::unordered_map<const RooAbsReal*,std::pair<int,int>> nServersClients;
    // used for dynamically scheduling each step's computations
    std::queue<const RooAbsReal*> computeQueue;
    std::unordered_map<const RooAbsReal*,std::pair<int,int>> nRemainingServersClients;
    std::queue<double*> buffers;
};

#endif //ROO_FIT_DRIVER_H
