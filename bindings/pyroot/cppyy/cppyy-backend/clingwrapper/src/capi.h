#ifndef CPPYY_CAPI
#define CPPYY_CAPI

#include <stddef.h>
#include <stdint.h>
#include "precommondefs.h"

#ifdef __cplusplus
extern "C" {
#endif // ifdef __cplusplus

    typedef size_t        cppyy_scope_t;
    typedef cppyy_scope_t cppyy_type_t;
    typedef void*         cppyy_object_t;
    typedef intptr_t      cppyy_method_t;

    typedef size_t        cppyy_index_t;
    typedef void*         cppyy_funcaddr_t;

    typedef unsigned long cppyy_exctype_t;

    /* direct interpreter access ---------------------------------------------- */
    RPY_EXPORTED
    int cppyy_compile(const char* code);

    /* name to opaque C++ scope representation -------------------------------- */
    RPY_EXPORTED
    char* cppyy_resolve_name(const char* cppitem_name);
    RPY_EXPORTED
    char* cppyy_resolve_enum(const char* enum_type);
    RPY_EXPORTED
    cppyy_scope_t cppyy_get_scope(const char* scope_name);
    RPY_EXPORTED
    cppyy_type_t cppyy_actual_class(cppyy_type_t klass, cppyy_object_t obj);
    RPY_EXPORTED
    size_t cppyy_size_of_klass(cppyy_type_t klass);
    RPY_EXPORTED
    size_t cppyy_size_of_type(const char* type_name);

    /* memory management ------------------------------------------------------ */
    RPY_EXPORTED
    cppyy_object_t cppyy_allocate(cppyy_type_t type);
    RPY_EXPORTED
    void cppyy_deallocate(cppyy_type_t type, cppyy_object_t self);
    RPY_EXPORTED
    cppyy_object_t cppyy_construct(cppyy_type_t type);
    RPY_EXPORTED
    void cppyy_destruct(cppyy_type_t type, cppyy_object_t self);

    /* method/function dispatching -------------------------------------------- */
    RPY_EXPORTED
    void cppyy_call_v(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);
    RPY_EXPORTED
    unsigned char cppyy_call_b(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);
    RPY_EXPORTED
    char cppyy_call_c(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);
    RPY_EXPORTED
    short cppyy_call_h(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);
    RPY_EXPORTED
    int cppyy_call_i(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);
    RPY_EXPORTED
    long cppyy_call_l(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);
    RPY_EXPORTED
    long long cppyy_call_ll(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);
    RPY_EXPORTED
    float cppyy_call_f(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);
    RPY_EXPORTED
    double cppyy_call_d(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);
    RPY_EXPORTED
    long double cppyy_call_ld(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);
    RPY_EXPORTED
    double cppyy_call_nld(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);

    RPY_EXPORTED
    void* cppyy_call_r(cppyy_method_t method, cppyy_object_t self, int nargs, void* args);
    RPY_EXPORTED
    char* cppyy_call_s(cppyy_method_t method, cppyy_object_t self, int nargs, void* args, size_t* length);
    RPY_EXPORTED
    cppyy_object_t cppyy_constructor(cppyy_method_t method, cppyy_type_t klass, int nargs, void* args);
    RPY_EXPORTED
    void cppyy_destructor(cppyy_type_t type, cppyy_object_t self);
    RPY_EXPORTED
    cppyy_object_t cppyy_call_o(cppyy_method_t method, cppyy_object_t self, int nargs, void* args, cppyy_type_t result_type);

    RPY_EXPORTED
    cppyy_funcaddr_t cppyy_function_address(cppyy_method_t method);

    /* handling of function argument buffer ----------------------------------- */
    RPY_EXPORTED
    void* cppyy_allocate_function_args(int nargs);
    RPY_EXPORTED
    void cppyy_deallocate_function_args(void* args);
    RPY_EXPORTED
    size_t cppyy_function_arg_sizeof();
    RPY_EXPORTED
    size_t cppyy_function_arg_typeoffset();

    /* scope reflection information ------------------------------------------- */
    RPY_EXPORTED
    int cppyy_is_namespace(cppyy_scope_t scope);
    RPY_EXPORTED
    int cppyy_is_template(const char* template_name);
    RPY_EXPORTED
    int cppyy_is_abstract(cppyy_type_t type);
    RPY_EXPORTED
    int cppyy_is_enum(const char* type_name);

    RPY_EXPORTED
    const char** cppyy_get_all_cpp_names(cppyy_scope_t scope, size_t* count);

    /* namespace reflection information --------------------------------------- */
    RPY_EXPORTED
    cppyy_index_t* cppyy_get_using_namespaces(cppyy_scope_t scope);

    /* class reflection information ------------------------------------------- */
    RPY_EXPORTED
    char* cppyy_final_name(cppyy_type_t type);
    RPY_EXPORTED
    char* cppyy_scoped_final_name(cppyy_type_t type);
    RPY_EXPORTED
    int cppyy_has_virtual_destructor(cppyy_type_t type);
    RPY_EXPORTED
    int cppyy_has_complex_hierarchy(cppyy_type_t type);
    RPY_EXPORTED
    int cppyy_num_bases(cppyy_type_t type);
    RPY_EXPORTED
    char* cppyy_base_name(cppyy_type_t type, int base_index);
    RPY_EXPORTED
    int cppyy_is_subtype(cppyy_type_t derived, cppyy_type_t base);
    RPY_EXPORTED
    int cppyy_is_smartptr(cppyy_type_t type);
    RPY_EXPORTED
    int cppyy_smartptr_info(const char* name, cppyy_type_t* raw, cppyy_method_t* deref);
    RPY_EXPORTED
    void cppyy_add_smartptr_type(const char* type_name);

    /* calculate offsets between declared and actual type, up-cast: direction > 0; down-cast: direction < 0 */
    RPY_EXPORTED
    ptrdiff_t cppyy_base_offset(cppyy_type_t derived, cppyy_type_t base, cppyy_object_t address, int direction);

    /* method/function reflection information --------------------------------- */
    RPY_EXPORTED
    int cppyy_num_methods(cppyy_scope_t scope);
    RPY_EXPORTED
    cppyy_index_t* cppyy_method_indices_from_name(cppyy_scope_t scope, const char* name);

    RPY_EXPORTED
    cppyy_method_t cppyy_get_method(cppyy_scope_t scope, cppyy_index_t idx);

    RPY_EXPORTED
    char* cppyy_method_name(cppyy_method_t);
    RPY_EXPORTED
    char* cppyy_method_full_name(cppyy_method_t);
    RPY_EXPORTED
    char* cppyy_method_mangled_name(cppyy_method_t);
    RPY_EXPORTED
    char* cppyy_method_result_type(cppyy_method_t);
    RPY_EXPORTED
    int cppyy_method_num_args(cppyy_method_t);
    RPY_EXPORTED
    int cppyy_method_req_args(cppyy_method_t);
    RPY_EXPORTED
    char* cppyy_method_arg_name(cppyy_method_t,int arg_index);
    RPY_EXPORTED
    char* cppyy_method_arg_type(cppyy_method_t, int arg_index);
    RPY_EXPORTED
    char* cppyy_method_arg_default(cppyy_method_t, int arg_index);
    RPY_EXPORTED
    char* cppyy_method_signature(cppyy_method_t, int show_formalargs);
    RPY_EXPORTED
    char* cppyy_method_signature_max(cppyy_method_t, int show_formalargs, int maxargs);
    RPY_EXPORTED
    char* cppyy_method_prototype(cppyy_scope_t scope, cppyy_method_t, int show_formalargs);
    RPY_EXPORTED
    int cppyy_is_const_method(cppyy_method_t);

    RPY_EXPORTED
    int cppyy_get_num_templated_methods(cppyy_scope_t scope);
    RPY_EXPORTED
    char* cppyy_get_templated_method_name(cppyy_scope_t scope, cppyy_index_t imeth);
    RPY_EXPORTED
    int cppyy_exists_method_template(cppyy_scope_t scope, const char* name);
    RPY_EXPORTED
    int cppyy_method_is_template(cppyy_scope_t scope, cppyy_index_t idx);
    RPY_EXPORTED
    cppyy_method_t cppyy_get_method_template(cppyy_scope_t scope, const char* name, const char* proto);

    RPY_EXPORTED
    cppyy_index_t cppyy_get_global_operator(
        cppyy_scope_t scope, cppyy_scope_t lc, cppyy_scope_t rc, const char* op);

    /* method properties ------------------------------------------------------ */
    RPY_EXPORTED
    int cppyy_is_publicmethod(cppyy_method_t);
    RPY_EXPORTED
    int cppyy_is_constructor(cppyy_method_t);
    RPY_EXPORTED
    int cppyy_is_destructor(cppyy_method_t);
    RPY_EXPORTED
    int cppyy_is_staticmethod(cppyy_method_t);

    /* data member reflection information ------------------------------------- */
    RPY_EXPORTED
    int cppyy_num_datamembers(cppyy_scope_t scope);
    RPY_EXPORTED
    char* cppyy_datamember_name(cppyy_scope_t scope, int datamember_index);
    RPY_EXPORTED
    char* cppyy_datamember_type(cppyy_scope_t scope, int datamember_index);
    RPY_EXPORTED
    intptr_t cppyy_datamember_offset(cppyy_scope_t scope, int datamember_index);
    RPY_EXPORTED
    int cppyy_datamember_index(cppyy_scope_t scope, const char* name);

    /* data member properties ------------------------------------------------- */
    RPY_EXPORTED
    int cppyy_is_publicdata(cppyy_type_t type, cppyy_index_t datamember_index);
    RPY_EXPORTED
    int cppyy_is_staticdata(cppyy_type_t type, cppyy_index_t datamember_index);
    RPY_EXPORTED
    int cppyy_is_const_data(cppyy_scope_t scope, cppyy_index_t idata);
    RPY_EXPORTED
    int cppyy_is_enum_data(cppyy_scope_t scope, cppyy_index_t idata);
    RPY_EXPORTED
    int cppyy_get_dimension_size(cppyy_scope_t scope, cppyy_index_t idata, int dimension);

    /* misc helpers ----------------------------------------------------------- */
    RPY_EXPORTED
    long long cppyy_strtoll(const char* str);
    RPY_EXPORTED
    unsigned long long cppyy_strtoull(const char* str);
    RPY_EXPORTED
    void cppyy_free(void* ptr);

    RPY_EXPORTED
    cppyy_object_t cppyy_charp2stdstring(const char* str, size_t sz);
    RPY_EXPORTED
    const char* cppyy_stdstring2charp(cppyy_object_t ptr, size_t* lsz);
    RPY_EXPORTED
    cppyy_object_t cppyy_stdstring2stdstring(cppyy_object_t ptr);

    RPY_EXPORTED
    double cppyy_longdouble2double(void*);
    RPY_EXPORTED
    void   cppyy_double2longdouble(double, void*);

    RPY_EXPORTED
    int         cppyy_vectorbool_getitem(cppyy_object_t ptr, int idx);
    RPY_EXPORTED
    void        cppyy_vectorbool_setitem(cppyy_object_t ptr, int idx, int value);

#ifdef __cplusplus
}
#endif // ifdef __cplusplus

#endif // ifndef CPPYY_CAPI
