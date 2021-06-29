// @(#)root/tmva/pymva $Id$
// Author: Sanjiban Sengupta, 2021

#include "TMVA/PyConverters.h"


#include "TMVA/Types.h"
#include "Rtypes.h"
#include "TString.h"


#include "TMVA/RModel.hxx"
#include "TMVA/RModelParser_Keras.h"


namespace TMVA{
namespace Experimental{
namespace SOFIE{

namespace PyKeras {
void ConvertToRoot(std::string filepath)
{
   RModel rmodel= SOFIE::PyKeras::Parse(filepath);
}
}

}
}
}

