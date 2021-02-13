#ifndef STD_STREAMS_H 
#define STD_STREAMS_H 1

#ifndef __CINT__
#include <ios>
#endif
#include <iostream>

#ifndef __CINT__
extern template class std::basic_ios<char,std::char_traits<char> >;
#endif

#endif // STD_STREAMS_H
