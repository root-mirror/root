#include "templates.h"


// template methods
long MyTemplatedMethodClass::get_size() { return -1; }

long MyTemplatedMethodClass::get_char_size()   { return (long)sizeof(char); }
long MyTemplatedMethodClass::get_int_size()    { return (long)sizeof(int); }
long MyTemplatedMethodClass::get_long_size()   { return (long)42; /* "lying" */ }
long MyTemplatedMethodClass::get_float_size()  { return (long)sizeof(float); }
long MyTemplatedMethodClass::get_double_size() { return (long)sizeof(double); }
long MyTemplatedMethodClass::get_self_size()   { return (long)sizeof(MyTemplatedMethodClass); }


// variadic templates
#ifdef WIN32
__declspec(dllexport)
#endif
std::string some_variadic::gTypeName = "";


// template with empty body
namespace T_WithEmptyBody {

#ifdef WIN32
__declspec(dllexport)
#endif
std::string side_effect = "not set";

template<typename T>
void some_empty() {
    side_effect = "side effect";
}

template void some_empty<int>();

} // namespace T_WithRValue
