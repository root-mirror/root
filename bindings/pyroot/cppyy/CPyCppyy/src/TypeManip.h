#ifndef CPYCPPYY_TYPEMANIP_H
#define CPYCPPYY_TYPEMANIP_H

#include <string>
#include <vector>


namespace CPyCppyy {

namespace TypeManip {

    std::string remove_const(const std::string& cppname);
    std::string clean_type(const std::string& cppname,
            bool template_strip = true, bool const_strip = true);
    std::string template_base(const std::string& cppname);

    void cppscope_to_pyscope(std::string& cppscope);
    std::string extract_namespace(const std::string& name);

    std::vector<std::string> extract_arg_types(const std::string& sig);

} // namespace TypeManip

} // namespace CPyCppyy

#endif // !CPYCPPYY_TYPEMANIP_H
