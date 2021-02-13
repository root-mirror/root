%module swig_runvector
%{
#include "runvector.h"
%}

%include "runvector.h"
%include "std_vector.i"

namespace std {
   %template(vectori) vector<int>;
}

