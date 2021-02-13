#ifndef ROOT_strtok
#define ROOT_strtok

#include <ROOT/RConfig.hxx>

#include <cstring>

// On Windows strtok_r does not exist: the function is called strtok_s.
inline char *R__STRTOK_R(char *str, const char *delim, char **saveptr)
{
#if defined(R__WIN32)
   return strtok_s(str, delim, saveptr);
#else
   return strtok_r(str, delim, saveptr);
#endif
}

#endif