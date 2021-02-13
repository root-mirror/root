#include "TError.h"

int main()
{
#ifdef NDEBUG
   Error("checkAssertsNDEBUG", "Compiling without assertions (NDEBUG flag is enabled) but asserts are configured");
   return 1;
#else
   return 0;
#endif
}
