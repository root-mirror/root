// define the TRandomGen classes
#ifdef DEFINE_TEMPL_INSTANCE

#include "TRandomGen.h"
#include "Math/MixMaxEngine.h"
#include "Math/RanluxppEngine.h"
#include "Math/StdEngine.h"

// define the instance
class  TRandomGen<ROOT::Math::MixMaxEngine<240,0>>;
class TRandomGen<ROOT::Math::MixMaxEngine<256,2>>; 
class TRandomGen<ROOT::Math::MixMaxEngine<256,4>>; 
class TRandomGen<ROOT::Math::MixMaxEngine<17,0>>;
class TRandomGen<ROOT::Math::MixMaxEngine<17,1>>;

class TRandomGen<ROOT::Math::RanluxppEngine2048>;

class  TRandomGen<ROOT::Math::StdEngine<std::mt19937_64> >;
class  TRandomGen<ROOT::Math::StdEngine<std::ranlux48> >;

#endif
