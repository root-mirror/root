// Bindings
#include "CPyCppyy.h"
#include "TypeManip.h"

// Standard
#include <ctype.h>


//- helpers ------------------------------------------------------------------
static inline
bool is_varchar(char c) {
    return isalnum((int)c) || c == '_' || c == ')' || c == '(' /* for (anonymous) */;
}

static inline
std::string::size_type find_qualifier_index(const std::string& name)
{
// Find the first location that is not part of the class name proper.
    std::string::size_type i = name.size() - 1;
    for ( ; 0 < i; --i) {
        std::string::value_type c = name[i];
        if (is_varchar(c) || c == '>') {
            if (c == 't' && 6 < i && !is_varchar(name[i-5]) && name.substr(i-4, 5) == "const")
                i -= 4;      // this skips 'const' on a pointer type
            else
                break;
        }
    }

    return i+1;
}

static inline void erase_const(std::string& name)
{
// Find and remove all occurrence of 'const'.
    std::string::size_type spos = std::string::npos;
    std::string::size_type start = 0;
    while ((spos = name.find("const", start)) != std::string::npos) {
    // make sure not to erase 'const' as part of the name: if it is
    // connected, before or after, to a variable name, then keep it
        std::string::size_type after = spos+5;
        if (after < name.size() && is_varchar(name[after])) {
            start = after;
            continue;
        } else if (after == name.size()) {
            if (spos > 0 && is_varchar(name[spos - 1]))
                break;
        }

        std::string::size_type i = 5;
        while (name[spos+i] == ' ') ++i;
        name.swap(name.erase(spos, i));
    }
}

static inline void rstrip(std::string& name)
{
// Remove space from the right side of name.
    std::string::size_type i = name.size();
    for ( ; 0 < i; --i) {
       if (!isspace(name[i]))
           break;
    }

    if (i != name.size())
        name = name.substr(0, i);
}


//----------------------------------------------------------------------------
std::string CPyCppyy::TypeManip::remove_const(const std::string& cppname)
{
// Remove 'const' qualifiers from the given C++ name.
    std::string::size_type tmplt_start = cppname.find('<');
    std::string::size_type type_stop   = cppname.rfind('>');
    if (cppname.find("::", type_stop+1) != std::string::npos) // e.g. klass<T>::some_typedef
        type_stop = cppname.find(' ', type_stop+1);
    if (tmplt_start != std::string::npos) {
    // only replace const qualifying cppname, not in template parameters
        std::string pre = cppname.substr(0, tmplt_start);
        erase_const(pre);
        std::string post = "";
        if (type_stop != std::string::npos) {
            post = cppname.substr(type_stop+1, std::string::npos);
            erase_const(post);
        }

        type_stop = type_stop == std::string::npos ? std::string::npos : type_stop+1;
        return pre + cppname.substr(tmplt_start, type_stop) + post;
    }

    std::string clean_name = cppname;
    erase_const(clean_name);
    return clean_name;
}

//----------------------------------------------------------------------------
std::string CPyCppyy::TypeManip::clean_type(
    const std::string& cppname, bool template_strip, bool const_strip)
{
// Strip C++ name from all qualifiers and compounds.
    std::string::size_type i = find_qualifier_index(cppname);
    std::string name = cppname.substr(0, i);
    rstrip(name);

    if (name.back() == ']') {                      // array type?
    // TODO: this fails templates instantiated on arrays (not common)
        name = name.substr(0, name.find('['));
    } else if (template_strip && name.back() == '>') {
        name = name.substr(0, name.find('<'));
    }

    if (const_strip) {
        if (template_strip)
            erase_const(name);
        else
            name = remove_const(name);
    }

    return name;
}

//----------------------------------------------------------------------------
std::string CPyCppyy::TypeManip::template_base(const std::string& cppname)
{
// If this is a template, return the underlying template name w/o arguments
    if (cppname.empty() || cppname.back() != '>')
        return cppname;

    int tpl_open = 0;
    for (std::string::size_type pos = cppname.size()-1; 0 < pos; --pos) {
        std::string::value_type c = cppname[pos];

    // count '<' and '>' to be able to skip template contents
        if (c == '>')
            ++tpl_open;
        else if (c == '<')
            --tpl_open;

        if (tpl_open == 0)
            return cppname.substr(0, pos);
    }

    return cppname;
}

//----------------------------------------------------------------------------
void CPyCppyy::TypeManip::cppscope_to_pyscope(std::string& cppscope)
{
// Change '::' in C++ scope into '.' as in a Python scope.
    std::string::size_type pos = 0;
    while ((pos = cppscope.find("::", pos)) != std::string::npos) {
        cppscope.replace(pos, 2, ".");
        pos += 1;
    }
}

//----------------------------------------------------------------------------
std::string CPyCppyy::TypeManip::extract_namespace(const std::string& name)
{
// Find the namespace the named class lives in, take care of templates
    if (name.empty())
        return name;

    int tpl_open = 0;
    for (std::string::size_type pos = name.size()-1; 0 < pos; --pos) {
        std::string::value_type c = name[pos];

    // count '<' and '>' to be able to skip template contents
        if (c == '>')
            ++tpl_open;
        else if (c == '<')
            --tpl_open;

    // collect name up to "::"
        else if (tpl_open == 0 && c == ':' && name[pos-1] == ':') {
        // found the extend of the scope ... done
            return name.substr(0, pos-1);
        }
    }

// no namespace; assume outer scope
    return "";
}

//----------------------------------------------------------------------------
std::vector<std::string> CPyCppyy::TypeManip::extract_arg_types(const std::string& sig)
{
// break out the argument types from the signature string
    std::vector<std::string> result;

    if (sig.empty() || sig == "()")
        return result;

    int tpl_open = 0;
    std::string::size_type start = 1;
    for (std::string::size_type pos = 1; pos < sig.size()-1; ++pos) {
        std::string::value_type c = sig[pos];

    // count '<' and '>' to be able to skip template contents
        if (c == '>')
            ++tpl_open;
        else if (c == '<')
            --tpl_open;

    // collect type name up to ',' or end ')'
        else if (tpl_open == 0 && c == ',') {
        // found the extend of the scope ... done
            result.push_back(sig.substr(start, pos-start));
            start = pos+1;
        }
    }

// add last type
    result.push_back(sig.substr(start, sig.rfind(")")-start));

    return result;
}

