#ifndef ROO_FIT_DRIVER_H
#define ROO_FIT_DRIVER_H

#include "rbc.h"

#include <queue>
#include <unordered_map>

class RooAbsData;
class RooAbsArg;
class RooAbsReal;
class RooArgSet;
class RooNLLVarNew;

class RooFitDriver {
  public:
     RooFitDriver(const RooAbsData& data, const RooNLLVarNew& topNode, int batchMode);
     ~RooFitDriver();
     //~  inline RooAbsReal* getTopNode() { return initialQueue.back(); }
     double getVal();
     RooArgSet* getParameters();
     
    struct NodeInfo {
        int nServers = 0;
        int nClients = 0;
    };

  private:

    double* getAvailableBuffer();

    const int batchMode=0;
    double* cudaMemDataset=nullptr;
    // used for preserving static info about the computation graph
    rbc::DataMap dataMap;
    const RooNLLVarNew& topNode;
    size_t nEvents;

    RooAbsData const* _data = nullptr;

    std::queue<const RooAbsReal*> initialQueue;

    std::unordered_map<const RooAbsArg*,NodeInfo> _nodeInfos;
    // used for dynamically scheduling each step's computations
    std::queue<const RooAbsReal*> computeQueue;

    std::queue<double*> _vectorBuffers;
};

#endif //ROO_FIT_DRIVER_H
