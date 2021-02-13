#include "Rtypes.h"

class TwoFloats {
public:
   virtual ~TwoFloats() {} // to make dictionary generation happy
   float a = 0;
   float b = 0;
   TwoFloats(){};
   TwoFloats(float d) : a(d), b(2 * d){};
   void Set(float d)
   {
      a = d;
      b = 2 * d;
   }
   ClassDef(TwoFloats, 2)
};

