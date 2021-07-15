// @(#)root/tmva/pymva $Id$
// Author: Sanjiban Sengupta, 2021

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Functionality for parsing a saved PyTorch .PT model into RModel object    *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Sanjiban Sengupta <sanjiban.sg@gmail.com> - IIIT, Bhubaneswar             *
 *                                                                                *
 * Copyright (c) 2020:                                                            *
 *      CERN, Switzerland                                                         *
 *      IIIT, Bhubaneswar                                                         *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/


#ifndef TMVA_SOFIE_RMODELPARSER_PYTORCH
#define TMVA_SOFIE_RMODELPARSER_PYTORCH

#include "TMVA/RModelParser_Common.h"

namespace TMVA{
namespace Experimental{
namespace SOFIE{

enum class NodeType{
   GEMM = 0, RELU = 1, TRANSPOSE = 2 //order sensitive

};

namespace PyTorch{
    RModel Parse(std::string filepath,std::vector<std::vector<size_t>> inputShapes, std::vector<ETensorType> dtype);
    RModel Parse(std::string filepath,std::vector<std::vector<size_t>> inputShapes);
  }
}
}
}

#endif //TMVA_PYMVA_RMODELPARSER_PYTORCH
