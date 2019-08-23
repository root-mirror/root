#ifndef __JIT_FUNCTIONS_HXX_
#define __JIT_FUNCTIONS_HXX_
//#include "unique_bdt.h"

#include "TInterpreter.h" // for gInterpreter
#include "bdt_helpers.hxx"

//////////////////////////////////////////////////////
/// JITTING FUNCTIONS
//////////////////////////////////////////////////////

/// JIT forest from sringed code
std::function<bool(const float *event)> jit_branchless_forest(std::string tojit, const std::string s_namespace = "")
{
   gInterpreter->Declare(tojit.c_str());
   bool use_namespaces = (!s_namespace.empty());

   std::string func_ref_name;
   if (use_namespaces)
      func_ref_name = "#pragma cling optimize(3)\n & branchless_" + s_namespace + "::branchless_generated_forest";
   else
      func_ref_name = "#pragma cling optimize(3)\n & branchless_generated_forest";
   auto ptr                    = gInterpreter->Calc(func_ref_name.c_str());
   bool (*func)(const float *) = reinterpret_cast<bool (*)(const float *)>(ptr);
   std::function<bool(const float *)> fWrapped{func};
   return fWrapped;
}

/// JIT forest from sringed code
template <typename T>
std::function<bool(const std::vector<float> &)> jit_branched_forest(std::string tojit, std::string s_namespace = "")
{
   gInterpreter->Declare(tojit.c_str());
   bool use_namespaces = (!s_namespace.empty());

   std::string func_ref_name;
   if (use_namespaces)
      // func_ref_name = "&s_f_" + s_namespace+"::generated_forest";
      func_ref_name = "#pragma cling optimize(3)\n &s_f_" + s_namespace + "::generated_forest";
   else
      func_ref_name = "#pragma cling optimize(3)\n &generated_forest";
   auto ptr                                 = gInterpreter->Calc(func_ref_name.c_str());
   bool (*func)(const std::vector<float> &) = reinterpret_cast<bool (*)(const std::vector<float> &)>(ptr);
   std::function<bool(const std::vector<float> &)> fWrapped{func};
   return fWrapped;
}

#endif
