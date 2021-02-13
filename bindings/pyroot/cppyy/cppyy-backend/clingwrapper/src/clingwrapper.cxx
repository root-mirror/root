// Bindings
#include "capi.h"
#include "cpp_cppyy.h"
#include "callcontext.h"

// ROOT
#include "TBaseClass.h"
#include "TClass.h"
#include "TClassRef.h"
#include "TClassTable.h"
#include "TClassEdit.h"
#include "TCollection.h"
#include "TDataMember.h"
#include "TDataType.h"
#include "TEnum.h"
#include "TEnumConstant.h"
#include "TEnv.h"
#include "TError.h"
#include "TException.h"
#include "TFunction.h"
#include "TFunctionTemplate.h"
#include "TGlobal.h"
#include "THashList.h"
#include "TInterpreter.h"
#include "TList.h"
#include "TListOfDataMembers.h"
#include "TListOfEnums.h"
#include "TMethod.h"
#include "TMethodArg.h"
#include "TROOT.h"
#include "TSystem.h"

// Standard
#include <assert.h>
#include <algorithm>     // for std::count, std::remove
#include <stdexcept>
#include <map>
#include <new>
#include <set>
#include <sstream>
#include <signal.h>
#include <stdlib.h>      // for getenv
#include <string.h>
#include <typeinfo>

// temp
#include <iostream>
typedef CPyCppyy::Parameter Parameter;
// --temp


// small number that allows use of stack for argument passing
const int SMALL_ARGS_N = 8;

// data for life time management ---------------------------------------------
typedef std::vector<TClassRef> ClassRefs_t;
static ClassRefs_t g_classrefs(1);
static const ClassRefs_t::size_type GLOBAL_HANDLE = 1;
static const ClassRefs_t::size_type STD_HANDLE = GLOBAL_HANDLE + 1;

typedef std::map<std::string, ClassRefs_t::size_type> Name2ClassRefIndex_t;
static Name2ClassRefIndex_t g_name2classrefidx;

namespace {

static inline Cppyy::TCppType_t find_memoized(const std::string& name)
{
    auto icr = g_name2classrefidx.find(name);
    if (icr != g_name2classrefidx.end())
        return (Cppyy::TCppType_t)icr->second;
    return (Cppyy::TCppType_t)0;
}

class CallWrapper {
public:
    typedef const void* DeclId_t;

public:
    CallWrapper(TFunction* f) : fDecl(f->GetDeclId()), fName(f->GetName()), fTF(nullptr) {}
    CallWrapper(DeclId_t fid, const std::string& n) : fDecl(fid), fName(n), fTF(nullptr) {}
    ~CallWrapper() {
        if (fTF && fDecl == fTF->GetDeclId())
            delete fTF;
    }

public:
    TInterpreter::CallFuncIFacePtr_t   fFaceptr;
    DeclId_t      fDecl;
    std::string   fName;
    TFunction*    fTF;
};

}

static std::vector<CallWrapper*> gWrapperHolder;
static inline CallWrapper* new_CallWrapper(TFunction* f)
{
    CallWrapper* wrap = new CallWrapper(f);
    gWrapperHolder.push_back(wrap);
    return wrap;
}

static inline CallWrapper* new_CallWrapper(CallWrapper::DeclId_t fid, const std::string& n)
{
    CallWrapper* wrap = new CallWrapper(fid, n);
    gWrapperHolder.push_back(wrap);
    return wrap;
}

typedef std::vector<TGlobal*> GlobalVars_t;
static GlobalVars_t g_globalvars;

static std::set<std::string> gSTLNames;


// data ----------------------------------------------------------------------
Cppyy::TCppScope_t Cppyy::gGlobalScope = GLOBAL_HANDLE;

// builtin types (including a few common STL templates as long as they live in
// the global namespace b/c of choices upstream)
static std::set<std::string> g_builtins =
    {"bool", "char", "signed char", "unsigned char", "wchar_t", "short", "unsigned short",
     "int", "unsigned int", "long", "unsigned long", "long long", "unsigned long long",
     "float", "double", "long double", "void",
     "allocator", "array", "basic_string", "complex", "initializer_list", "less", "list",
     "map", "pair", "set", "vector"};

// smart pointer types
static std::set<std::string> gSmartPtrTypes =
    {"auto_ptr", "std::auto_ptr", "shared_ptr", "std::shared_ptr",
     "unique_ptr", "std::unique_ptr", "weak_ptr", "std::weak_ptr"};

// to filter out ROOT names
static std::set<std::string> gInitialNames;
static std::set<std::string> gRootSOs;

// configuration
static bool gEnableFastPath = true;


// global initialization -----------------------------------------------------
namespace {

// names copied from TUnixSystem
#ifdef WIN32
const int SIGBUS   = 0;      // simple placeholders for ones that don't exist
const int SIGSYS   = 0;
const int SIGPIPE  = 0;
const int SIGQUIT  = 0;
const int SIGWINCH = 0;
const int SIGALRM  = 0;
const int SIGCHLD  = 0;
const int SIGURG   = 0;
const int SIGUSR1  = 0;
const int SIGUSR2  = 0;
#endif

static struct Signalmap_t {
   int               fCode;
   const char       *fSigName;
} gSignalMap[kMAXSIGNALS] = {       // the order of the signals should be identical
   { SIGBUS,   "bus error" }, // to the one in TSysEvtHandler.h
   { SIGSEGV,  "segmentation violation" },
   { SIGSYS,    "bad argument to system call" },
   { SIGPIPE,   "write on a pipe with no one to read it" },
   { SIGILL,    "illegal instruction" },
   { SIGABRT,   "abort" },
   { SIGQUIT,   "quit" },
   { SIGINT,    "interrupt" },
   { SIGWINCH,  "window size change" },
   { SIGALRM,   "alarm clock" },
   { SIGCHLD,   "death of a child" },
   { SIGURG,    "urgent data arrived on an I/O channel" },
   { SIGFPE,    "floating point exception" },
   { SIGTERM,   "termination signal" },
   { SIGUSR1,   "user-defined signal 1" },
   { SIGUSR2,   "user-defined signal 2" }
};

static void inline do_trace(int sig) {
    std::cerr << " *** Break *** " << (sig < kMAXSIGNALS ? gSignalMap[sig].fSigName : "") << std::endl;
    gSystem->StackTrace();
}

class TExceptionHandlerImp : public TExceptionHandler {
public:
    virtual void HandleException(Int_t sig) {
        if (TROOT::Initialized()) {
            if (gException) {
                gInterpreter->RewindDictionary();
                gInterpreter->ClearFileBusy();
            }

            if (!getenv("CPPYY_CRASH_QUIET"))
                do_trace(sig);

        // jump back, if catch point set
            Throw(sig);
        }

        do_trace(sig);
        gSystem->Exit(128 + sig);
    }
};

class ApplicationStarter {
public:
    ApplicationStarter() {
    // initialize ROOT early to guarantee proper order of shutdown later on (gROOT is a
    // macro that resolves to the ROOT::GetROOT() function call)
        (void)gROOT;

    // setup dummy holders for global and std namespaces
        assert(g_classrefs.size() == GLOBAL_HANDLE);
        g_name2classrefidx[""]      = GLOBAL_HANDLE;
        g_classrefs.push_back(TClassRef(""));

    // aliases for std (setup already in pythonify)
        g_name2classrefidx["std"]   = STD_HANDLE;
        g_name2classrefidx["::std"] = g_name2classrefidx["std"];
        g_classrefs.push_back(TClassRef("std"));

    // add a dummy global to refer to as null at index 0
        g_globalvars.push_back(nullptr);

    // disable fast path if requested
        if (getenv("CPPYY_DISABLE_FASTPATH")) gEnableFastPath = false;

    // fill the set of STL names
        const char* stl_names[] = {"allocator", "auto_ptr", "bad_alloc", "bad_cast",
            "bad_exception", "bad_typeid", "basic_filebuf", "basic_fstream", "basic_ifstream",
            "basic_ios", "basic_iostream", "basic_istream", "basic_istringstream",
            "basic_ofstream", "basic_ostream", "basic_ostringstream", "basic_streambuf",
            "basic_string", "basic_stringbuf", "basic_stringstream", "binary_function",
            "binary_negate", "bitset", "byte", "char_traits", "codecvt_byname", "codecvt", "collate",
            "collate_byname", "compare", "complex", "ctype_byname", "ctype", "default_delete",
            "deque", "divides", "domain_error", "equal_to", "exception", "forward_list", "fpos",
            "function", "greater_equal", "greater", "gslice_array", "gslice", "hash", "indirect_array",
            "integer_sequence", "invalid_argument", "ios_base", "istream_iterator", "istreambuf_iterator",
            "istrstream", "iterator_traits", "iterator", "length_error", "less_equal", "less",
            "list", "locale", "localedef utility", "locale utility", "logic_error", "logical_and",
            "logical_not", "logical_or", "map", "mask_array", "mem_fun", "mem_fun_ref", "messages",
            "messages_byname", "minus", "modulus", "money_get", "money_put", "moneypunct",
            "moneypunct_byname", "multimap", "multiplies", "multiset", "negate", "not_equal_to",
            "num_get", "num_put", "numeric_limits", "numpunct", "numpunct_byname",
            "ostream_iterator", "ostreambuf_iterator", "ostrstream", "out_of_range",
            "overflow_error", "pair", "plus", "pointer_to_binary_function",
            "pointer_to_unary_function", "priority_queue", "queue", "range_error",
            "raw_storage_iterator", "reverse_iterator", "runtime_error", "set", "shared_ptr",
            "slice_array", "slice", "stack", "string", "strstream", "strstreambuf",
            "time_get_byname", "time_get", "time_put_byname", "time_put", "unary_function",
            "unary_negate", "unique_ptr", "underflow_error", "unordered_map", "unordered_multimap",
            "unordered_multiset", "unordered_set", "valarray", "vector", "weak_ptr", "wstring"};
        for (auto& name : stl_names)
            gSTLNames.insert(name);

    // set opt level (default to 2 if not given; Cling itself defaults to 0)
        int optLevel = 2;
        if (getenv("CPPYY_OPT_LEVEL")) optLevel = atoi(getenv("CPPYY_OPT_LEVEL"));
        if (optLevel != 0) {
            std::ostringstream s;
            s << "#pragma cling optimize " << optLevel;
            gInterpreter->ProcessLine(s.str().c_str());
        }

    // load frequently used headers
        const char* code =
               "#include <iostream>\n"
               "#include <string>\n"
               "#include <DllImport.h>\n"     // defines R__EXTERN
               "#include <vector>\n"
               "#include <utility>";
        gInterpreter->ProcessLine(code);

    // make sure we run in batch mode as far as ROOT graphics is concerned
        if (!getenv("ROOTSYS"))
            gROOT->SetBatch(kTRUE);

    // create helpers for comparing thingies
        gInterpreter->Declare(
            "namespace __cppyy_internal { template<class C1, class C2>"
            " bool is_equal(const C1& c1, const C2& c2) { return (bool)(c1 == c2); } }");
        gInterpreter->Declare(
            "namespace __cppyy_internal { template<class C1, class C2>"
            " bool is_not_equal(const C1& c1, const C2& c2) { return (bool)(c1 != c2); } }");

    // retrieve all initial (ROOT) C++ names in the global scope to allow filtering later
        if (!getenv("CPPYY_NO_ROOT_FILTER")) {
            gROOT->GetListOfGlobals(true);             // force initialize
            gROOT->GetListOfGlobalFunctions(true);     // id.
            std::set<std::string> initial;
            Cppyy::GetAllCppNames(GLOBAL_HANDLE, initial);
            gInitialNames = initial;

#ifndef WIN32
            gRootSOs.insert("libCore.so ");
            gRootSOs.insert("libRIO.so ");
            gRootSOs.insert("libThread.so ");
            gRootSOs.insert("libMathCore.so ");
#else
            gRootSOs.insert("libCore.dll ");
            gRootSOs.insert("libRIO.dll ");
            gRootSOs.insert("libThread.dll ");
            gRootSOs.insert("libMathCore.dll ");
#endif
        }

    // start off with a reasonable size placeholder for wrappers
        gWrapperHolder.reserve(1024);

    // create an exception handler to process signals
        gExceptionHandler = new TExceptionHandlerImp{};
    }

    ~ApplicationStarter() {
        for (auto wrap : gWrapperHolder)
            delete wrap;
        delete gExceptionHandler; gExceptionHandler = nullptr;
    }
} _applicationStarter;

} // unnamed namespace


// local helpers -------------------------------------------------------------
static inline
TClassRef& type_from_handle(Cppyy::TCppScope_t scope)
{
    assert((ClassRefs_t::size_type)scope < g_classrefs.size());
    return g_classrefs[(ClassRefs_t::size_type)scope];
}

static inline
TFunction* m2f(Cppyy::TCppMethod_t method) {
    CallWrapper* wrap = ((CallWrapper*)method);
    if (!wrap->fTF || wrap->fTF->GetDeclId() != wrap->fDecl) {
        MethodInfo_t* mi = gInterpreter->MethodInfo_Factory(wrap->fDecl);
        wrap->fTF = new TFunction(mi);
    }
    return wrap->fTF;
}

static inline
char* cppstring_to_cstring(const std::string& cppstr)
{
    char* cstr = (char*)malloc(cppstr.size()+1);
    memcpy(cstr, cppstr.c_str(), cppstr.size()+1);
    return cstr;
}

static inline
bool match_name(const std::string& tname, const std::string fname)
{
// either match exactly, or match the name as template
    if (fname.rfind(tname, 0) == 0) {
        if ((tname.size() == fname.size()) ||
              (tname.size() < fname.size() && fname[tname.size()] == '<'))
           return true;
    }
    return false;
}

static inline
bool is_missclassified_stl(const std::string& name)
{
    std::string::size_type pos = name.find('<');
    if (pos != std::string::npos)
        return gSTLNames.find(name.substr(0, pos)) != gSTLNames.end();
    return gSTLNames.find(name) != gSTLNames.end();
}


// direct interpreter access -------------------------------------------------
bool Cppyy::Compile(const std::string& code)
{
    return gInterpreter->Declare(code.c_str());
}


// name to opaque C++ scope representation -----------------------------------
std::string Cppyy::ResolveName(const std::string& cppitem_name)
{
// Fully resolve the given name to the final type name.

// try memoized type cache, in case seen before
    TCppType_t klass = find_memoized(cppitem_name);
    if (klass) return GetScopedFinalName(klass);

// remove global scope '::' if present
    std::string tclean = cppitem_name.compare(0, 2, "::") == 0 ?
        cppitem_name.substr(2, std::string::npos) : cppitem_name;

// classes (most common)
    tclean = TClassEdit::CleanType(tclean.c_str());
    if (tclean.empty() /* unknown, eg. an operator */) return cppitem_name;

// reduce [N] to []
    if (tclean[tclean.size()-1] == ']')
        tclean = tclean.substr(0, tclean.rfind('[')) + "[]";

    if (tclean.rfind("byte", 0) == 0 || tclean.rfind("std::byte", 0) == 0)
        return tclean;

// check data types list (accept only builtins as typedefs will
// otherwise not be resolved)
    TDataType* dt = gROOT->GetType(tclean.c_str());
    if (dt && dt->GetType() != kOther_t) return dt->GetFullTypeName();

// special case for enums
    if (IsEnum(cppitem_name))
        return ResolveEnum(cppitem_name);

// special case for clang's builtin __type_pack_element (which does not resolve)
    if (cppitem_name.rfind("__type_pack_element", 0) != std::string::npos) {
    // shape is "__type_pack_element<index,type1,type2,...,typeN>cpd": extract
    // first the index, and from there the indexed type; finally, restore the
    // qualifiers
        const char* str = cppitem_name.c_str();
        char* endptr = nullptr;
        unsigned long index = strtoul(str+20, &endptr, 0);

        std::string tmplvars{endptr};
        auto start = tmplvars.find(',') + 1;
        auto end = tmplvars.find(',', start);
        while (index != 0) {
            start = end+1;
            end = tmplvars.find(',', start);
            if (end == std::string::npos) end = tmplvars.rfind('>');
            --index;
        }

        std::string resolved = tmplvars.substr(start, end-start);
        auto cpd = tmplvars.rfind('>');
        if (cpd != std::string::npos && cpd+1 != tmplvars.size())
            return resolved + tmplvars.substr(cpd+1, std::string::npos);
        return resolved;
    }

// typedefs
    return TClassEdit::ResolveTypedef(tclean.c_str(), true);
}

static std::map<std::string, std::string> resolved_enum_types;
std::string Cppyy::ResolveEnum(const std::string& enum_type)
{
// The underlying type of a an enum may be any kind of integer.
// Resolve that type via a workaround (note: this function assumes
// that the enum_type name is a valid enum type name)
    auto res = resolved_enum_types.find(enum_type);
    if (res != resolved_enum_types.end())
        return res->second;

// desugar the type before resolving
    std::string et_short = TClassEdit::ShortType(enum_type.c_str(), 1);
    if (et_short.find("(anonymous") == std::string::npos) {
        std::ostringstream decl;
    // TODO: now presumed fixed with https://sft.its.cern.ch/jira/browse/ROOT-6988
        for (auto& itype : {"unsigned int"}) {
            decl << "std::is_same<"
                 << itype
                 << ", std::underlying_type<"
                 << et_short
                 << ">::type>::value;";
            if (gInterpreter->ProcessLine(decl.str().c_str())) {
            // TODO: "re-sugaring" like this is brittle, but the top
            // should be re-translated into AST-based code anyway
                std::string resugared;
                if (et_short.size() != enum_type.size()) {
                    auto pos = enum_type.find(et_short);
                    if (pos != std::string::npos) {
                        resugared = enum_type.substr(0, pos) + itype;
                        if (pos+et_short.size() < enum_type.size())
                            resugared += enum_type.substr(pos+et_short.size(), std::string::npos);
                    }
                }
                if (resugared.empty()) resugared = itype;
                resolved_enum_types[enum_type] = resugared;
                return resugared;
            }
        }
    }

// failed or anonymous ... signal upstream to special case this
    int ipos = (int)enum_type.size()-1;
    for (; 0 <= ipos; --ipos) {
        char c = enum_type[ipos];
        if (isspace(c)) continue;
        if (isalnum(c) || c == '_' || c == '>' || c == ')') break;
    }
    bool isConst = enum_type.find("const ", 6) != std::string::npos;
    std::string restype = isConst ? "const " : "";
    restype += "internal_enum_type_t"+enum_type.substr((std::string::size_type)ipos+1, std::string::npos);
    resolved_enum_types[enum_type] = restype;
    return restype;     // should default to some int variant
}

Cppyy::TCppScope_t Cppyy::GetScope(const std::string& sname)
{
// First, try cache
    TCppType_t result = find_memoized(sname);
    if (result) return result;

// Second, skip builtins before going through the more expensive steps of resolving
// typedefs and looking up TClass
    if (g_builtins.find(sname) != g_builtins.end())
        return (TCppScope_t)0;

// TODO: scope_name should always be final already?
// Resolve name fully before lookup to make sure all aliases point to the same scope
    std::string scope_name = ResolveName(sname);
    bool bHasAlias = sname != scope_name;
    if (bHasAlias) {
        result = find_memoized(scope_name);
        if (result) return result;
    }

// both failed, but may be STL name that's missing 'std::' now, but didn't before
    bool b_scope_name_missclassified = is_missclassified_stl(scope_name);
    if (b_scope_name_missclassified) {
        result = find_memoized("std::"+scope_name);
        if (result) g_name2classrefidx["std::"+scope_name] = (ClassRefs_t::size_type)result;
    }
    bool b_sname_missclassified = bHasAlias ? is_missclassified_stl(sname) : false;
    if (b_sname_missclassified) {
        if (!result) result = find_memoized("std::"+sname);
        if (result) g_name2classrefidx["std::"+sname] = (ClassRefs_t::size_type)result;
    }

    if (result) return result;

// use TClass directly, to enable auto-loading; class may be stubbed (eg. for
// function returns) or forward declared, leading to a non-null TClass that is
// otherwise invalid/unusable
    TClassRef cr(TClass::GetClass(scope_name.c_str(), true /* load */, true /* silent */));
    if (!cr.GetClass())
        return (TCppScope_t)0;

// memoize found/created TClass
    ClassRefs_t::size_type sz = g_classrefs.size();
    g_name2classrefidx[scope_name] = sz;
    if (bHasAlias) g_name2classrefidx[sname] = sz;
    g_classrefs.push_back(TClassRef(scope_name.c_str()));

// TODO: make ROOT/meta NOT remove std :/
    if (b_scope_name_missclassified)
        g_name2classrefidx["std::"+scope_name] = sz;
    if (b_sname_missclassified)
        g_name2classrefidx["std::"+sname] = sz;

    return (TCppScope_t)sz;
}

bool Cppyy::IsTemplate(const std::string& template_name)
{
    return (bool)gInterpreter->CheckClassTemplate(template_name.c_str());
}

namespace {
    class AutoCastRTTI {
    public:
        virtual ~AutoCastRTTI() {}
    };
}

Cppyy::TCppType_t Cppyy::GetActualClass(TCppType_t klass, TCppObject_t obj)
{
    TClassRef& cr = type_from_handle(klass);
    if (!cr.GetClass() || !obj) return klass;

#ifdef _WIN64
// Cling does not provide a consistent ImageBase address for calculating relative addresses
// as used in Windows 64b RTTI. So, check for our own RTTI extension instead. If that fails,
// see whether the unmangled raw_name is available (e.g. if this is an MSVC compiled rather
// than JITed class) and pass on if it is.
    volatile const char* raw = nullptr;     // to prevent too aggressive reordering
    try {
    // this will filter those objects that do not have RTTI to begin with (throws)
        AutoCastRTTI* pcst = (AutoCastRTTI*)obj;
        raw = typeid(*pcst).raw_name();

    // check the signature id (0 == absolute, 1 == relative, 2 == ours)
        void* vfptr = *(void**)((intptr_t)obj);
        void* meta = (void*)((intptr_t)*((void**)((intptr_t)vfptr-sizeof(void*))));
        if (*(intptr_t*)meta == 2) {
        // access the extra data item which is an absolute pointer to the RTTI
            void* ptdescr = (void*)((intptr_t)meta + 4*sizeof(unsigned long)+sizeof(void*));
            if (ptdescr && *(void**)ptdescr) {
                auto rtti = *(std::type_info**)ptdescr;
                raw = rtti->raw_name();
                if (raw && raw[0] != '\0')    // likely unnecessary
                    return (TCppType_t)GetScope(rtti->name());
            }

            return klass;        // do not fall through if no RTTI info available
        }

    // if the raw name is the empty string (no guarantees that this is so as truly, the
    // address is corrupt, but it is common to be empty), then there is no accessible RTTI
    // and getting the unmangled name will crash ...
        if (!raw || raw[0] == '\0')
            return klass;
    } catch (std::bad_typeid) {
        return klass;        // can't risk passing to ROOT/meta as it may do RTTI
    }
#endif

    TClass* clActual = cr->GetActualClass((void*)obj);
    if (clActual && clActual != cr.GetClass()) {
        auto itt = g_name2classrefidx.find(clActual->GetName());
        if (itt != g_name2classrefidx.end())
            return (TCppType_t)itt->second;
        return (TCppType_t)GetScope(clActual->GetName());
    }

    return klass;
}

size_t Cppyy::SizeOf(TCppType_t klass)
{
    TClassRef& cr = type_from_handle(klass);
    if (cr.GetClass() && cr->GetClassInfo())
        return (size_t)gInterpreter->ClassInfo_Size(cr->GetClassInfo());
    return (size_t)0;
}

size_t Cppyy::SizeOf(const std::string& type_name)
{
    TDataType* dt = gROOT->GetType(type_name.c_str());
    if (dt) return dt->Size();
    return SizeOf(GetScope(type_name));
}

bool Cppyy::IsBuiltin(const std::string& type_name)
{
    TDataType* dt = gROOT->GetType(TClassEdit::CleanType(type_name.c_str(), 1).c_str());
    if (dt) return dt->GetType() != kOther_t;
    return false;
}

bool Cppyy::IsComplete(const std::string& type_name)
{
// verify whether the dictionary of this class is fully available
    bool b = false;

    int oldEIL = gErrorIgnoreLevel;
    gErrorIgnoreLevel = 3000;
    TClass* klass = TClass::GetClass(TClassEdit::ShortType(type_name.c_str(), 1).c_str());
    if (klass && klass->GetClassInfo())     // works for normal case w/ dict
        b = gInterpreter->ClassInfo_IsLoaded(klass->GetClassInfo());
    else {    // special case for forward declared classes
        ClassInfo_t* ci = gInterpreter->ClassInfo_Factory(type_name.c_str());
        if (ci) {
            b = gInterpreter->ClassInfo_IsLoaded(ci);
            gInterpreter->ClassInfo_Delete(ci);  // we own the fresh class info
        }
    }
    gErrorIgnoreLevel = oldEIL;
    return b;
}

// memory management ---------------------------------------------------------
Cppyy::TCppObject_t Cppyy::Allocate(TCppType_t type)
{
    TClassRef& cr = type_from_handle(type);
    return (TCppObject_t)malloc(gInterpreter->ClassInfo_Size(cr->GetClassInfo()));
}

void Cppyy::Deallocate(TCppType_t /* type */, TCppObject_t instance)
{
    ::operator delete(instance);
}

Cppyy::TCppObject_t Cppyy::Construct(TCppType_t type)
{
    TClassRef& cr = type_from_handle(type);
    return (TCppObject_t)cr->New();
}

static std::map<Cppyy::TCppType_t, bool> sHasOperatorDelete;
void Cppyy::Destruct(TCppType_t type, TCppObject_t instance)
{
    TClassRef& cr = type_from_handle(type);
    if (cr->ClassProperty() & (kClassHasExplicitDtor | kClassHasImplicitDtor))
        cr->Destructor((void*)instance);
    else {
        ROOT::DelFunc_t fdel = cr->GetDelete();
        if (fdel) fdel((void*)instance);
        else {
            auto ib = sHasOperatorDelete.find(type);
            if (ib == sHasOperatorDelete.end()) {
                sHasOperatorDelete[type] = (bool)cr->GetListOfAllPublicMethods()->FindObject("operator delete");
                ib = sHasOperatorDelete.find(type);
            }
            ib->second ? cr->Destructor((void*)instance) : free((void*)instance);
        }
    }
}


// method/function dispatching -----------------------------------------------
static TInterpreter::CallFuncIFacePtr_t GetCallFunc(Cppyy::TCppMethod_t method)
{
// TODO: method should be a callfunc, so that no mapping would be needed.
    CallWrapper* wrap = (CallWrapper*)method;

    CallFunc_t* callf = gInterpreter->CallFunc_Factory();
    MethodInfo_t* meth = gInterpreter->MethodInfo_Factory(wrap->fDecl);
    gInterpreter->CallFunc_SetFunc(callf, meth);
    gInterpreter->MethodInfo_Delete(meth);

    if (!(callf && gInterpreter->CallFunc_IsValid(callf))) {
    // TODO: propagate this error to caller w/o use of Python C-API
    /*
        PyErr_Format(PyExc_RuntimeError, "could not resolve %s::%s(%s)",
            const_cast<TClassRef&>(klass).GetClassName(),
            wrap.fName, callString.c_str()); */
        std::cerr << "TODO: report unresolved function error to Python\n";
        if (callf) gInterpreter->CallFunc_Delete(callf);
        return TInterpreter::CallFuncIFacePtr_t{};
    }

// generate the wrapper and JIT it; ignore wrapper generation errors (will simply
// result in a nullptr that is reported upstream if necessary; often, however,
// there is a different overload available that will do)
    auto oldErrLvl = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;
    wrap->fFaceptr = gInterpreter->CallFunc_IFacePtr(callf);
    gErrorIgnoreLevel = oldErrLvl;

    gInterpreter->CallFunc_Delete(callf);   // does not touch IFacePtr
    return wrap->fFaceptr;
}

static inline
bool copy_args(Parameter* args, size_t nargs, void** vargs)
{
    bool runRelease = false;
    for (size_t i = 0; i < nargs; ++i) {
        switch (args[i].fTypeCode) {
        case 'X':       /* (void*)type& with free */
            runRelease = true;
        case 'V':       /* (void*)type& */
            vargs[i] = args[i].fValue.fVoidp;
            break;
        case 'r':       /* const type& */
            vargs[i] = args[i].fRef;
            break;
        default:        /* all other types in union */
            vargs[i] = (void*)&args[i].fValue.fVoidp;
            break;
        }
    }
    return runRelease;
}

static inline
void release_args(Parameter* args, size_t nargs) {
    for (size_t i = 0; i < nargs; ++i) {
        if (args[i].fTypeCode == 'X')
            free(args[i].fValue.fVoidp);
    }
}

static inline bool WrapperCall(Cppyy::TCppMethod_t method, size_t nargs, void* args_, void* self, void* result)
{
    Parameter* args = (Parameter*)args_;

    CallWrapper* wrap = (CallWrapper*)method;
    const TInterpreter::CallFuncIFacePtr_t& faceptr = wrap->fFaceptr.fGeneric ? wrap->fFaceptr : GetCallFunc(method);
    if (!faceptr.fGeneric)
        return false;        // happens with compilation error

    if (faceptr.fKind == TInterpreter::CallFuncIFacePtr_t::kGeneric) {
        bool runRelease = false;
        if (nargs <= SMALL_ARGS_N) {
            void* smallbuf[SMALL_ARGS_N];
            if (nargs) runRelease = copy_args(args, nargs, smallbuf);
            faceptr.fGeneric(self, (int)nargs, smallbuf, result);
        } else {
            std::vector<void*> buf(nargs);
            runRelease = copy_args(args, nargs, buf.data());
            faceptr.fGeneric(self, (int)nargs, buf.data(), result);
        }
        if (runRelease) release_args(args, nargs);
        return true;
    }

    if (faceptr.fKind == TInterpreter::CallFuncIFacePtr_t::kCtor) {
        bool runRelease = false;
        if (nargs <= SMALL_ARGS_N) {
            void* smallbuf[SMALL_ARGS_N];
            if (nargs) runRelease = copy_args(args, nargs, (void**)smallbuf);
            faceptr.fCtor((void**)smallbuf, result, (unsigned long)nargs);
        } else {
            std::vector<void*> buf(nargs);
            runRelease = copy_args(args, nargs, buf.data());
            faceptr.fCtor(buf.data(), result, (unsigned long)nargs);
        }
        if (runRelease) release_args(args, nargs);
        return true;
    }

    if (faceptr.fKind == TInterpreter::CallFuncIFacePtr_t::kDtor) {
        std::cerr << " DESTRUCTOR NOT IMPLEMENTED YET! " << std::endl;
        return false;
    }

    return false;
}

template<typename T>
static inline
T CallT(Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, size_t nargs, void* args)
{
    T t{};
    if (WrapperCall(method, nargs, args, (void*)self, &t))
        return t;
    return (T)-1;
}

#define CPPYY_IMP_CALL(typecode, rtype)                                      \
rtype Cppyy::Call##typecode(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args)\
{                                                                            \
    return CallT<rtype>(method, self, nargs, args);                          \
}

void Cppyy::CallV(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args)
{
    if (!WrapperCall(method, nargs, args, (void*)self, nullptr))
        return /* TODO ... report error */;
}

CPPYY_IMP_CALL(B,  unsigned char)
CPPYY_IMP_CALL(C,  char         )
CPPYY_IMP_CALL(H,  short        )
CPPYY_IMP_CALL(I,  int          )
CPPYY_IMP_CALL(L,  long         )
CPPYY_IMP_CALL(LL, Long64_t     )
CPPYY_IMP_CALL(F,  float        )
CPPYY_IMP_CALL(D,  double       )
CPPYY_IMP_CALL(LD, LongDouble_t )

void* Cppyy::CallR(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args)
{
    void* r = nullptr;
    if (WrapperCall(method, nargs, args, (void*)self, &r))
        return r;
    return nullptr;
}

char* Cppyy::CallS(
    TCppMethod_t method, TCppObject_t self, size_t nargs, void* args, size_t* length)
{
    char* cstr = nullptr;
    TClassRef cr("std::string");
    std::string* cppresult = (std::string*)malloc(sizeof(std::string));
    if (WrapperCall(method, nargs, args, self, (void*)cppresult)) {
        cstr = cppstring_to_cstring(*cppresult);
        *length = cppresult->size();
        cppresult->std::string::~basic_string();
    } else
        *length = 0;
    free((void*)cppresult);
    return cstr;
}

Cppyy::TCppObject_t Cppyy::CallConstructor(
    TCppMethod_t method, TCppType_t /* klass */, size_t nargs, void* args)
{
    void* obj = nullptr;
    if (WrapperCall(method, nargs, args, nullptr, &obj))
        return (TCppObject_t)obj;
    return (TCppObject_t)0;
}

void Cppyy::CallDestructor(TCppType_t type, TCppObject_t self)
{
    TClassRef& cr = type_from_handle(type);
    cr->Destructor((void*)self, true);
}

Cppyy::TCppObject_t Cppyy::CallO(TCppMethod_t method,
    TCppObject_t self, size_t nargs, void* args, TCppType_t result_type)
{
    TClassRef& cr = type_from_handle(result_type);
    void* obj = ::operator new(gInterpreter->ClassInfo_Size(cr->GetClassInfo()));
    if (WrapperCall(method, nargs, args, self, obj))
        return (TCppObject_t)obj;
    return (TCppObject_t)0;
}

Cppyy::TCppFuncAddr_t Cppyy::GetFunctionAddress(TCppMethod_t method, bool check_enabled)
{
    if (check_enabled && !gEnableFastPath) return (TCppFuncAddr_t)nullptr;
    TFunction* f = m2f(method);
    return (TCppFuncAddr_t)gInterpreter->FindSym(f->GetMangledName());
}


// handling of function argument buffer --------------------------------------
void* Cppyy::AllocateFunctionArgs(size_t nargs)
{
    return new Parameter[nargs];
}

void Cppyy::DeallocateFunctionArgs(void* args)
{
    delete [] (Parameter*)args;
}

size_t Cppyy::GetFunctionArgSizeof()
{
    return sizeof(Parameter);
}

size_t Cppyy::GetFunctionArgTypeoffset()
{
    return offsetof(Parameter, fTypeCode);
}


// scope reflection information ----------------------------------------------
bool Cppyy::IsNamespace(TCppScope_t scope)
{
// Test if this scope represents a namespace.
    if (scope == GLOBAL_HANDLE)
        return true;
    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass())
        return cr->Property() & kIsNamespace;
    return false;
}

bool Cppyy::IsAbstract(TCppType_t klass)
{
// Test if this type may not be instantiated.
    TClassRef& cr = type_from_handle(klass);
    if (cr.GetClass())
        return cr->Property() & kIsAbstract;
    return false;
}

bool Cppyy::IsEnum(const std::string& type_name)
{
    if (type_name.empty()) return false;
    std::string tn_short = TClassEdit::ShortType(type_name.c_str(), 1);
    if (tn_short.empty()) return false;
    return gInterpreter->ClassInfo_IsEnum(tn_short.c_str());
}

// helpers for stripping scope names
static
std::string outer_with_template(const std::string& name)
{
// Cut down to the outer-most scope from <name>, taking proper care of templates.
    int tpl_open = 0;
    for (std::string::size_type pos = 0; pos < name.size(); ++pos) {
        std::string::value_type c = name[pos];

    // count '<' and '>' to be able to skip template contents
        if (c == '<')
            ++tpl_open;
        else if (c == '>')
            --tpl_open;

    // collect name up to "::"
        else if (tpl_open == 0 && \
                 c == ':' && pos+1 < name.size() && name[pos+1] == ':') {
        // found the extend of the scope ... done
            return name.substr(0, pos-1);
        }
    }

// whole name is apparently a single scope
    return name;
}

static
std::string outer_no_template(const std::string& name)
{
// Cut down to the outer-most scope from <name>, drop templates
    std::string::size_type first_scope = name.find(':');
    if (first_scope == std::string::npos)
        return name.substr(0, name.find('<'));
    std::string::size_type first_templ = name.find('<');
    if (first_templ == std::string::npos)
        return name.substr(0, first_scope);
    return name.substr(0, std::min(first_templ, first_scope));
}

#define FILL_COLL(type, filter) {                                             \
    TIter itr{coll};                                                          \
    type* obj = nullptr;                                                      \
    while ((obj = (type*)itr.Next())) {                                       \
        const char* nm = obj->GetName();                                      \
        if (nm && nm[0] != '_' && !(obj->Property() & (filter))) {            \
            if (gInitialNames.find(nm) == gInitialNames.end())                \
                cppnames.insert(nm);                                          \
    }}}

static inline
void cond_add(Cppyy::TCppScope_t scope, const std::string& ns_scope,
    std::set<std::string>& cppnames, const char* name, bool nofilter = false)
{
    if (!name || name[0] == '_' || strstr(name, ".h") != 0 || strncmp(name, "operator", 8) == 0)
        return;

    if (scope == GLOBAL_HANDLE) {
        std::string to_add = outer_no_template(name);
        if ((nofilter || gInitialNames.find(to_add) == gInitialNames.end()) && !is_missclassified_stl(name))
            cppnames.insert(outer_no_template(name));
    } else if (scope == STD_HANDLE) {
        if (strncmp(name, "std::", 5) == 0) {
            name += 5;
#ifdef __APPLE__
            if (strncmp(name, "__1::", 5) == 0) name += 5;
#endif
        } else if (!is_missclassified_stl(name))
            return;
        cppnames.insert(outer_no_template(name));
    } else {
        if (strncmp(name, ns_scope.c_str(), ns_scope.size()) == 0)
            cppnames.insert(outer_with_template(name + ns_scope.size()));
    }
}

void Cppyy::GetAllCppNames(TCppScope_t scope, std::set<std::string>& cppnames)
{
// Collect all known names of C++ entities under scope. This is useful for IDEs
// employing tab-completion, for example. Note that functions names need not be
// unique as they can be overloaded.
    TClassRef& cr = type_from_handle(scope);
    if (scope != GLOBAL_HANDLE && !(cr.GetClass() && cr->Property()))
        return;

    std::string ns_scope = GetFinalName(scope);
    if (scope != GLOBAL_HANDLE) ns_scope += "::";

// add existing values from read rootmap files if within this scope
    TCollection* coll = gInterpreter->GetMapfile()->GetTable();
    {
        TIter itr{coll};
        TEnvRec* ev = nullptr;
        while ((ev = (TEnvRec*)itr.Next())) {
        // TEnv contains rootmap entries and user-side rootmap files may be already
        // loaded on startup. Thus, filter on file name rather than load time.
            if (gRootSOs.find(ev->GetValue()) == gRootSOs.end())
                cond_add(scope, ns_scope, cppnames, ev->GetName(), true);
        }
    }

// do we care about the class table or are the rootmap and list of types enough?
/*
    gClassTable->Init();
    const int N = gClassTable->Classes();
    for (int i = 0; i < N; ++i)
        cond_add(scope, ns_scope, cppnames, gClassTable->Next());
*/

// any other types (e.g. that may have come from parsing headers)
    coll = gROOT->GetListOfTypes();
    {
        TIter itr{coll};
        TDataType* dt = nullptr;
        while ((dt = (TDataType*)itr.Next())) {
            if (!(dt->Property() & kIsFundamental)) {
                cond_add(scope, ns_scope, cppnames, dt->GetName());
            }
        }
    }

// add functions
    coll = (scope == GLOBAL_HANDLE) ?
        gROOT->GetListOfGlobalFunctions() : cr->GetListOfMethods();
    {
        TIter itr{coll};
        TFunction* obj = nullptr;
        while ((obj = (TFunction*)itr.Next())) {
            const char* nm = obj->GetName();
        // skip templated functions, adding only the un-instantiated ones
            if (nm && nm[0] != '_' && strstr(nm, "<") == 0 && strncmp(nm, "operator", 8) != 0) {
                if (gInitialNames.find(nm) == gInitialNames.end())
                    cppnames.insert(nm);
            }
        }
    }

// add uninstantiated templates
    coll = (scope == GLOBAL_HANDLE) ?
        gROOT->GetListOfFunctionTemplates() : cr->GetListOfFunctionTemplates();
    FILL_COLL(TFunctionTemplate, kIsPrivate | kIsProtected)

// add (global) data members
    if (scope == GLOBAL_HANDLE) {
        coll = gROOT->GetListOfGlobals();
        FILL_COLL(TGlobal, kIsEnum | kIsPrivate | kIsProtected)
    } else {
        coll = cr->GetListOfDataMembers();
        FILL_COLL(TDataMember, kIsEnum | kIsPrivate | kIsProtected)
        coll = cr->GetListOfUsingDataMembers();
        FILL_COLL(TDataMember, kIsEnum | kIsPrivate | kIsProtected)
    }

// add enums values only for user classes/namespaces
    if (scope != GLOBAL_HANDLE && scope != STD_HANDLE) {
        coll = cr->GetListOfEnums();
        FILL_COLL(TEnum, kIsPrivate | kIsProtected)
    }

#ifdef __APPLE__
// special case for Apple, add version namespace '__1' entries to std
    if (scope == STD_HANDLE)
        GetAllCppNames(GetScope("std::__1"), cppnames);
#endif
}


// class reflection information ----------------------------------------------
std::vector<Cppyy::TCppScope_t> Cppyy::GetUsingNamespaces(TCppScope_t scope)
{
    std::vector<Cppyy::TCppScope_t> res;
    if (!IsNamespace(scope))
        return res;

#ifdef __APPLE__
    if (scope == STD_HANDLE) {
        res.push_back(GetScope("__1"));
        return res;
    }
#endif

    TClassRef& cr = type_from_handle(scope);
    if (!cr.GetClass() || !cr->GetClassInfo())
        return res;

    const std::vector<std::string>& v = gInterpreter->GetUsingNamespaces(cr->GetClassInfo());
    res.reserve(v.size());
    for (auto uid : v) {
        Cppyy::TCppScope_t uscope = GetScope(uid);
        if (uscope) res.push_back(uscope);
    }

    return res;
}


// class reflection information ----------------------------------------------
std::string Cppyy::GetFinalName(TCppType_t klass)
{
    if (klass == GLOBAL_HANDLE)
        return "";
    TClassRef& cr = type_from_handle(klass);
    std::string clName = cr->GetName();
// TODO: why is this template splitting needed?
    std::string::size_type pos = clName.substr(0, clName.find('<')).rfind("::");
    if (pos != std::string::npos)
        return clName.substr(pos+2, std::string::npos);
    return clName;
}

std::string Cppyy::GetScopedFinalName(TCppType_t klass)
{
    if (klass == GLOBAL_HANDLE)
        return "";
    TClassRef& cr = type_from_handle(klass);
    if (cr.GetClass()) {
        std::string name = cr->GetName();
        if (is_missclassified_stl(name))
            return std::string("std::")+cr->GetName();
        return cr->GetName();
    }
    return "";
}

bool Cppyy::HasVirtualDestructor(TCppType_t klass)
{
    TClassRef& cr = type_from_handle(klass);
    if (!cr.GetClass())
        return false;

    TFunction* f = cr->GetMethod(("~"+GetFinalName(klass)).c_str(), "");
    if (f && (f->Property() & kIsVirtual))
        return true;

    return false;
}

bool Cppyy::HasComplexHierarchy(TCppType_t klass)
{
    int is_complex = 1;
    size_t nbases = 0;

    TClassRef& cr = type_from_handle(klass);
    if (cr.GetClass() && cr->GetListOfBases() != 0)
        nbases = GetNumBases(klass);

    if (1 < nbases)
        is_complex = 1;
    else if (nbases == 0)
        is_complex = 0;
    else {         // one base class only
        TBaseClass* base = (TBaseClass*)cr->GetListOfBases()->At(0);
        if (base->Property() & kIsVirtualBase)
            is_complex = 1;       // TODO: verify; can be complex, need not be.
        else
            is_complex = HasComplexHierarchy(GetScope(base->GetName()));
    }

    return is_complex;
}

Cppyy::TCppIndex_t Cppyy::GetNumBases(TCppType_t klass)
{
// Get the total number of base classes that this class has.
    TClassRef& cr = type_from_handle(klass);
    if (cr.GetClass() && cr->GetListOfBases() != 0)
        return (TCppIndex_t)cr->GetListOfBases()->GetSize();
    return (TCppIndex_t)0;
}

std::string Cppyy::GetBaseName(TCppType_t klass, TCppIndex_t ibase)
{
    TClassRef& cr = type_from_handle(klass);
    return ((TBaseClass*)cr->GetListOfBases()->At((int)ibase))->GetName();
}

bool Cppyy::IsSubtype(TCppType_t derived, TCppType_t base)
{
    if (derived == base)
        return true;
    TClassRef& derived_type = type_from_handle(derived);
    TClassRef& base_type = type_from_handle(base);
    return derived_type->GetBaseClass(base_type) != 0;
}

bool Cppyy::IsSmartPtr(TCppType_t klass)
{
    TClassRef& cr = type_from_handle(klass);
    const std::string& tn = cr->GetName();
    if (gSmartPtrTypes.find(tn.substr(0, tn.find("<"))) != gSmartPtrTypes.end())
        return true;
    return false;
}

bool Cppyy::GetSmartPtrInfo(
    const std::string& tname, TCppType_t* raw, TCppMethod_t* deref)
{
    const std::string& rn = ResolveName(tname);
    if (gSmartPtrTypes.find(rn.substr(0, rn.find("<"))) != gSmartPtrTypes.end()) {
        if (!raw && !deref) return true;

        TClassRef& cr = type_from_handle(GetScope(tname));
        if (cr.GetClass()) {
            TFunction* func = cr->GetMethod("operator->", "");
            if (!func) {
                gInterpreter->UpdateListOfMethods(cr.GetClass());
                func = cr->GetMethod("operator->", "");
            }
            if (func) {
               if (deref) *deref = (TCppMethod_t)new_CallWrapper(func);
               if (raw) *raw = GetScope(TClassEdit::ShortType(
                   func->GetReturnTypeNormalizedName().c_str(), 1));
               return (!deref || *deref) && (!raw || *raw);
            }
        }
    }

    return false;
}

void Cppyy::AddSmartPtrType(const std::string& type_name)
{
    gSmartPtrTypes.insert(ResolveName(type_name));
}


// type offsets --------------------------------------------------------------
ptrdiff_t Cppyy::GetBaseOffset(TCppType_t derived, TCppType_t base,
    TCppObject_t address, int direction, bool rerror)
{
// calculate offsets between declared and actual type, up-cast: direction > 0; down-cast: direction < 0
    if (derived == base || !(base && derived))
        return (ptrdiff_t)0;

    TClassRef& cd = type_from_handle(derived);
    TClassRef& cb = type_from_handle(base);

    if (!cd.GetClass() || !cb.GetClass())
        return (ptrdiff_t)0;

    ptrdiff_t offset = -1;
    if (!(cd->GetClassInfo() && cb->GetClassInfo())) {     // gInterpreter requirement
    // would like to warn, but can't quite determine error from intentional
    // hiding by developers, so only cover the case where we really should have
    // had a class info, but apparently don't:
        if (cd->IsLoaded()) {
        // warn to allow diagnostics
            std::ostringstream msg;
            msg << "failed offset calculation between " << cb->GetName() << " and " << cd->GetName();
            // TODO: propagate this warning to caller w/o use of Python C-API
            // PyErr_Warn(PyExc_RuntimeWarning, const_cast<char*>(msg.str().c_str()));
            std::cerr << "Warning: " << msg.str() << '\n';
        }

    // return -1 to signal caller NOT to apply offset
        return rerror ? (ptrdiff_t)offset : 0;
    }

    offset = gInterpreter->ClassInfo_GetBaseOffset(
        cd->GetClassInfo(), cb->GetClassInfo(), (void*)address, direction > 0);
    if (offset == -1)   // Cling error, treat silently
        return rerror ? (ptrdiff_t)offset : 0;

    return (ptrdiff_t)(direction < 0 ? -offset : offset);
}


// method/function reflection information ------------------------------------
Cppyy::TCppIndex_t Cppyy::GetNumMethods(TCppScope_t scope)
{
    if (IsNamespace(scope))
        return (TCppIndex_t)0;     // enforce lazy

    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass() && cr->GetListOfMethods(true)) {
        Cppyy::TCppIndex_t nMethods = (TCppIndex_t)cr->GetListOfMethods(false)->GetSize();
        if (nMethods == (TCppIndex_t)0) {
            std::string clName = GetScopedFinalName(scope);
            if (clName.find('<') != std::string::npos) {
            // chicken-and-egg problem: TClass does not know about methods until
            // instantiation, so force it
                if (clName.find("std::", 0, 5) == std::string::npos && \
                        is_missclassified_stl(clName)) {
                // TODO: this is too simplistic for template arguments missing std::
                    clName = "std::" + clName;
                }
                std::ostringstream stmt;
                stmt << "template class " << clName << ";";
                gInterpreter->Declare(stmt.str().c_str());

            // now reload the methods
                return (TCppIndex_t)cr->GetListOfMethods(true)->GetSize();
            }
        }
        return nMethods;
    }

    return (TCppIndex_t)0;         // unknown class?
}

std::vector<Cppyy::TCppIndex_t> Cppyy::GetMethodIndicesFromName(
    TCppScope_t scope, const std::string& name)
{
    std::vector<TCppIndex_t> indices;
    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass()) {
        gInterpreter->UpdateListOfMethods(cr.GetClass());
        int imeth = 0;
        TFunction* func = nullptr;
        TIter next(cr->GetListOfMethods());
        while ((func = (TFunction*)next())) {
            if (match_name(name, func->GetName())) {
                if (func->Property() & kIsPublic)
                    indices.push_back((TCppIndex_t)imeth);
            }
            ++imeth;
        }
    } else if (scope == GLOBAL_HANDLE) {
        TCollection* funcs = gROOT->GetListOfGlobalFunctions(true);

    // tickle deserialization
        if (!funcs->FindObject(name.c_str()))
            return indices;

        TFunction* func = nullptr;
        TIter ifunc(funcs);
        while ((func = (TFunction*)ifunc.Next())) {
            if (match_name(name, func->GetName()))
                indices.push_back((TCppIndex_t)new_CallWrapper(func));
        }
    }

    return indices;
}

Cppyy::TCppMethod_t Cppyy::GetMethod(TCppScope_t scope, TCppIndex_t idx)
{
    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass()) {
        TFunction* f = (TFunction*)cr->GetListOfMethods(false)->At((int)idx);
        if (f) return (Cppyy::TCppMethod_t)new_CallWrapper(f);
        return (Cppyy::TCppMethod_t)nullptr;
    }

    assert(klass == (Cppyy::TCppType_t)GLOBAL_HANDLE);
    return (Cppyy::TCppMethod_t)idx;
}

std::string Cppyy::GetMethodName(TCppMethod_t method)
{
    if (method) {
        const std::string& name = ((CallWrapper*)method)->fName;

        if (name.compare(0, 8, "operator") != 0)
        // strip template instantiation part, if any
            return name.substr(0, name.find('<'));
        return name;
    }
    return "<unknown>";
}

std::string Cppyy::GetMethodFullName(TCppMethod_t method)
{
    if (method) {
        std::string name = ((CallWrapper*)method)->fName;
        name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
        return name;
    }
    return "<unknown>";
}

std::string Cppyy::GetMethodMangledName(TCppMethod_t method)
{
    if (method)
        return m2f(method)->GetMangledName();
    return "<unknown>";
}

std::string Cppyy::GetMethodResultType(TCppMethod_t method)
{
    if (method) {
        TFunction* f = m2f(method);
        if (f->ExtraProperty() & kIsConstructor)
            return "constructor";
        std::string restype = f->GetReturnTypeName();
        // TODO: this is ugly, but we can't use GetReturnTypeName() for ostreams
        // and maybe others, whereas GetReturnTypeNormalizedName() has proven to
        // be save in all cases (Note: 'int8_t' covers 'int8_t' and 'uint8_t')
        if (restype.find("int8_t") != std::string::npos)
            return restype;
        restype = f->GetReturnTypeNormalizedName();
        if (restype == "(lambda)") {
            std::ostringstream s;
            // TODO: what if there are parameters to the lambda?
            s << "__cling_internal::FT<decltype("
              << GetMethodFullName(method) << "(";
            for (Cppyy::TCppIndex_t i = 0; i < Cppyy::GetMethodNumArgs(method); ++i) {
                if (i != 0) s << ", ";
                s << Cppyy::GetMethodArgType(method, i) << "{}";
            }
            s << "))>::F";
            TClass* cl = TClass::GetClass(s.str().c_str());
            if (cl) return cl->GetName();
            // TODO: signal some type of error (or should that be upstream?
        }
        return restype;
    }
    return "<unknown>";
}

Cppyy::TCppIndex_t Cppyy::GetMethodNumArgs(TCppMethod_t method)
{
    if (method)
        return m2f(method)->GetNargs();
    return 0;
}

Cppyy::TCppIndex_t Cppyy::GetMethodReqArgs(TCppMethod_t method)
{
    if (method) {
        TFunction* f = m2f(method);
        return (TCppIndex_t)(f->GetNargs() - f->GetNargsOpt());
    }
    return (TCppIndex_t)0;
}

std::string Cppyy::GetMethodArgName(TCppMethod_t method, TCppIndex_t iarg)
{
    if (method) {
        TFunction* f = m2f(method);
        TMethodArg* arg = (TMethodArg*)f->GetListOfMethodArgs()->At((int)iarg);
        return arg->GetName();
    }
    return "<unknown>";
}

std::string Cppyy::GetMethodArgType(TCppMethod_t method, TCppIndex_t iarg)
{
    if (method) {
        TFunction* f = m2f(method);
        TMethodArg* arg = (TMethodArg*)f->GetListOfMethodArgs()->At((int)iarg);
        return arg->GetTypeNormalizedName();
    }
    return "<unknown>";
}

std::string Cppyy::GetMethodArgDefault(TCppMethod_t method, TCppIndex_t iarg)
{
    if (method) {
        TFunction* f = m2f(method);
        TMethodArg* arg = (TMethodArg*)f->GetListOfMethodArgs()->At((int)iarg);
        const char* def = arg->GetDefault();
        if (def)
            return def;
    }

    return "";
}

std::string Cppyy::GetMethodSignature(TCppMethod_t method, bool show_formalargs, TCppIndex_t maxargs)
{
    TFunction* f = m2f(method);
    if (f) {
        std::ostringstream sig;
        sig << "(";
        int nArgs = f->GetNargs();
        if (maxargs != (TCppIndex_t)-1) nArgs = std::min(nArgs, (int)maxargs);
        for (int iarg = 0; iarg < nArgs; ++iarg) {
            TMethodArg* arg = (TMethodArg*)f->GetListOfMethodArgs()->At(iarg);
            sig << arg->GetFullTypeName();
            if (show_formalargs) {
                const char* argname = arg->GetName();
                if (argname && argname[0] != '\0') sig << " " << argname;
                const char* defvalue = arg->GetDefault();
                if (defvalue && defvalue[0] != '\0') sig << " = " << defvalue;
            }
            if (iarg != nArgs-1) sig << (show_formalargs ? ", " : ",");
        }
        sig << ")";
        return sig.str();
    }
    return "<unknown>";
}

std::string Cppyy::GetMethodPrototype(TCppScope_t scope, TCppMethod_t method, bool show_formalargs)
{
    std::string scName = GetScopedFinalName(scope);
    TFunction* f = m2f(method);
    if (f) {
        std::ostringstream sig;
        sig << f->GetReturnTypeName() << " "
            << scName << "::" << f->GetName();
        sig << GetMethodSignature(method, show_formalargs);
        return sig.str();
    }
    return "<unknown>";
}

bool Cppyy::IsConstMethod(TCppMethod_t method)
{
    if (method) {
        TFunction* f = m2f(method);
        return f->Property() & kIsConstMethod;
    }
    return false;
}

Cppyy::TCppIndex_t Cppyy::GetNumTemplatedMethods(TCppScope_t scope)
{
    if (scope == (TCppScope_t)GLOBAL_HANDLE) {
        TCollection* coll = gROOT->GetListOfFunctionTemplates();
        if (coll) return (TCppIndex_t)coll->GetSize();
    } else {
        TClassRef& cr = type_from_handle(scope);
        if (cr.GetClass()) {
            TCollection* coll = cr->GetListOfFunctionTemplates(true);
            if (coll) return (TCppIndex_t)coll->GetSize();
        }
    }

// failure ...
    return (TCppIndex_t)0;
}

std::string Cppyy::GetTemplatedMethodName(TCppScope_t scope, TCppIndex_t imeth)
{
    if (scope == (TCppScope_t)GLOBAL_HANDLE)
        return ((THashList*)gROOT->GetListOfFunctionTemplates())->At((int)imeth)->GetName();
    else {
        TClassRef& cr = type_from_handle(scope);
        if (cr.GetClass())
            return cr->GetListOfFunctionTemplates(false)->At((int)imeth)->GetName();
    }

// failure ...
    assert(!"should not be called unless GetNumTemplatedMethods() succeeded");
    return "";
}

bool Cppyy::IsTemplatedConstructor(TCppScope_t scope, TCppIndex_t imeth)
{
    if (scope == (TCppScope_t)GLOBAL_HANDLE)
        return false;

    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass()) {
        TFunctionTemplate* f = (TFunctionTemplate*)cr->GetListOfFunctionTemplates(false)->At((int)imeth);
        return f->ExtraProperty() & kIsConstructor;
    }

    return false;
}

bool Cppyy::ExistsMethodTemplate(TCppScope_t scope, const std::string& name)
{
    if (scope == (TCppScope_t)GLOBAL_HANDLE)
        return (bool)gROOT->GetFunctionTemplate(name.c_str());
    else {
        TClassRef& cr = type_from_handle(scope);
        if (cr.GetClass())
            return (bool)cr->GetFunctionTemplate(name.c_str());
    }

// failure ...
    return false;
}

bool Cppyy::IsMethodTemplate(TCppScope_t scope, TCppIndex_t idx)
{
    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass()) {
        TFunction* f = (TFunction*)cr->GetListOfMethods(false)->At((int)idx);
        if (f && strstr(f->GetName(), "<")) return true;
        return false;
    }

    assert(scope == (Cppyy::TCppType_t)GLOBAL_HANDLE);
    if (((CallWrapper*)idx)->fName.find('<') != std::string::npos) return true;
    return false;
}

// helpers for Cppyy::GetMethodTemplate()
static std::map<TDictionary::DeclId_t, CallWrapper*> gMethodTemplates;

Cppyy::TCppMethod_t Cppyy::GetMethodTemplate(
    TCppScope_t scope, const std::string& name, const std::string& proto)
{
// There is currently no clean way of extracting a templated method out of ROOT/meta
// for a variety of reasons, none of them fundamental. The game played below is to
// first get any pre-existing functions already managed by ROOT/meta, but if that fails,
// to do an explicit lookup that ignores the prototype (i.e. the full name should be
// enough), and finally to ignore the template arguments part of the name as this fails
// in cling if there are default parameters.
// It would be possible to get the prototype from the created functions and use that to
// do a new lookup, after which ROOT/meta will manage the function. However, neither
// TFunction::GetPrototype() nor TFunction::GetSignature() is of the proper form, so
// we'll/ manage the new TFunctions instead and will assume that they are cached on the
// calling side to prevent multiple creations.
    TFunction* func = nullptr; ClassInfo_t* cl = nullptr;
    if (scope == (cppyy_scope_t)GLOBAL_HANDLE) {
        func = gROOT->GetGlobalFunctionWithPrototype(name.c_str(), proto.c_str());
        if (func && name.back() == '>' && name != func->GetName())
            func = nullptr;  // happens if implicit conversion matches the overload
    } else {
        TClassRef& cr = type_from_handle(scope);
        if (cr.GetClass()) {
            func = cr->GetMethodWithPrototype(name.c_str(), proto.c_str());
            if (!func) {
                cl = cr->GetClassInfo();
            // try base classes to cover a common 'using' case (TODO: this is stupid and misses
            // out on base classes; fix that with improved access to Cling)
                TCppIndex_t nbases = GetNumBases(scope);
                for (TCppIndex_t i = 0; i < nbases; ++i) {
                    TClassRef& base = type_from_handle(GetScope(GetBaseName(scope, i)));
                    if (base.GetClass()) {
                        func = base->GetMethodWithPrototype(name.c_str(), proto.c_str());
                        if (func) break;
                    }
                }
            }
        }
    }

    if (!func && name.back() == '>' && (cl || scope == (cppyy_scope_t)GLOBAL_HANDLE)) {
    // try again, ignoring proto in case full name is complete template
        auto declid = gInterpreter->GetFunction(cl, name.c_str());
        if (declid) {
             auto existing = gMethodTemplates.find(declid);
             if (existing == gMethodTemplates.end()) {
                 auto cw = new_CallWrapper(declid, name);
                 existing = gMethodTemplates.insert(std::make_pair(declid, cw)).first;
             }
             return (TCppMethod_t)existing->second;
        }
    }

    if (func) {
    // make sure we didn't match a non-templated overload
        if (func->ExtraProperty() & kIsTemplateSpec)
            return (TCppMethod_t)new_CallWrapper(func);

    // disregard this non-templated method as it will be considered when appropriate
        return (TCppMethod_t)nullptr;
    }

// try again with template arguments removed from name, if applicable
    if (name.back() == '>') {
        auto pos = name.find('<');
        if (pos != std::string::npos) {
            TCppMethod_t cppmeth = GetMethodTemplate(scope, name.substr(0, pos), proto);
            if (cppmeth) {
            // allow if requested template names match up to the result
                const std::string& alt = GetMethodFullName(cppmeth);
                if (name.size() < alt.size() && alt.find('<') == pos) {
                    const std::string& partial = name.substr(pos, name.size()-1-pos);
                    if (strncmp(partial.c_str(), alt.substr(pos, alt.size()-1-pos).c_str(), partial.size()) == 0)
                        return cppmeth;
                }
            }
        }
    }

// failure ...
    return (TCppMethod_t)nullptr;
}

static inline
std::string type_remap(const std::string& n1, const std::string& n2)
{
// Operator lookups of (C++ string, Python str) should succeeded, for the combos of
// string/str, wstring/str, string/unicode and wstring/unicode; since C++ does not have a
// operator+(std::string, std::wstring), we'll have to look up the same type and rely on
// the converters in CPyCppyy/_cppyy.
    if (n1 == "str") {
        if (n2 == "std::basic_string<wchar_t,std::char_traits<wchar_t>,std::allocator<wchar_t> >")
            return n2;                      // match like for like
        return "std::string";               // probably best bet
    } else if (n1 == "float")
        return "double";                    // debatable, but probably intended
    return n1;
}

Cppyy::TCppIndex_t Cppyy::GetGlobalOperator(
    TCppType_t scope, const std::string& lc, const std::string& rc, const std::string& opname)
{
// Find a global operator function with a matching signature; prefer by-ref, but
// fall back on by-value if that fails.
    std::string lcname1 = TClassEdit::CleanType(lc.c_str());
    const std::string& rcname = rc.empty() ? rc : type_remap(TClassEdit::CleanType(rc.c_str()), lcname1);
    const std::string& lcname = type_remap(lcname1, rcname);

    std::string proto = lcname + "&" + (rc.empty() ? rc : (", " + rcname + "&"));
    if (scope == (cppyy_scope_t)GLOBAL_HANDLE) {
        TFunction* func = gROOT->GetGlobalFunctionWithPrototype(opname.c_str(), proto.c_str());
        if (func) return (TCppIndex_t)new_CallWrapper(func);
        proto = lcname + (rc.empty() ? rc : (", " + rcname));
        func = gROOT->GetGlobalFunctionWithPrototype(opname.c_str(), proto.c_str());
        if (func) return (TCppIndex_t)new_CallWrapper(func);
    } else {
        TClassRef& cr = type_from_handle(scope);
        if (cr.GetClass()) {
            TFunction* func = cr->GetMethodWithPrototype(opname.c_str(), proto.c_str());
            if (func) return (TCppIndex_t)cr->GetListOfMethods()->IndexOf(func);
            proto = lcname + (rc.empty() ? rc : (", " + rcname));
            func = cr->GetMethodWithPrototype(opname.c_str(), proto.c_str());
            if (func) return (TCppIndex_t)cr->GetListOfMethods()->IndexOf(func);
        }
    }

// failure ...
    return (TCppIndex_t)-1;
}

// method properties ---------------------------------------------------------
bool Cppyy::IsPublicMethod(TCppMethod_t method)
{
    if (method) {
        TFunction* f = m2f(method);
        return f->Property() & kIsPublic;
    }
    return false;
}

bool Cppyy::IsProtectedMethod(TCppMethod_t method)
{
    if (method) {
        TFunction* f = m2f(method);
        return f->Property() & kIsProtected;
    }
    return false;
}

bool Cppyy::IsConstructor(TCppMethod_t method)
{
    if (method) {
        TFunction* f = m2f(method);
        return f->ExtraProperty() & kIsConstructor;
    }
    return false;
}

bool Cppyy::IsDestructor(TCppMethod_t method)
{
    if (method) {
        TFunction* f = m2f(method);
        return f->ExtraProperty() & kIsDestructor;
    }
    return false;
}

bool Cppyy::IsStaticMethod(TCppMethod_t method)
{
    if (method) {
        TFunction* f = m2f(method);
        return f->Property() & kIsStatic;
    }
    return false;
}

// data member reflection information ----------------------------------------
Cppyy::TCppIndex_t Cppyy::GetNumDatamembers(TCppScope_t scope)
{
    if (IsNamespace(scope))
        return (TCppIndex_t)0;     // enforce lazy

    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass()) {
        Cppyy::TCppIndex_t sum = 0;
        if (cr->GetListOfDataMembers())
            sum = cr->GetListOfDataMembers()->GetSize();
        if (cr->GetListOfUsingDataMembers())
            sum += cr->GetListOfUsingDataMembers()->GetSize();
        return sum;
    }

    return (TCppIndex_t)0;         // unknown class?
}

static TDataMember *GetDataMemberByIndex(TClassRef cr, int idata)
{
    if (!cr.GetClass() || !cr->GetListOfDataMembers())
        return nullptr;

    int numDMs = cr->GetListOfDataMembers()->GetSize();
    if ((int)idata < numDMs)
        return (TDataMember*)cr->GetListOfDataMembers()->At((int)idata);
    return (TDataMember*)cr->GetListOfUsingDataMembers()->At((int)idata - numDMs);
}

std::string Cppyy::GetDatamemberName(TCppScope_t scope, TCppIndex_t idata)
{
    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass()) {
        TDataMember *m = GetDataMemberByIndex(cr, (int)idata);
        return m->GetName();
    }
    assert(scope == GLOBAL_HANDLE);
    TGlobal* gbl = g_globalvars[idata];
    return gbl->GetName();
}

std::string Cppyy::GetDatamemberType(TCppScope_t scope, TCppIndex_t idata)
{
    if (scope == GLOBAL_HANDLE) {
        TGlobal* gbl = g_globalvars[idata];
        std::string fullType = gbl->GetFullTypeName();

        if ((int)gbl->GetArrayDim() > 1)
            fullType.append("*");
        else if ((int)gbl->GetArrayDim() == 1) {
            std::ostringstream s;
            s << '[' << gbl->GetMaxIndex(0) << ']' << std::ends;
            fullType.append(s.str());
        }
        return fullType;
    }

    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass())  {
        TDataMember* m = GetDataMemberByIndex(cr, (int)idata);
    // TODO: fix this upstream. Usually, we want m->GetFullTypeName(), because it does
    // not resolve typedefs, but it looses scopes for inner classes/structs, so in that
    // case m->GetTrueTypeName() should be used (this also cleans up the cases where
    // the "full type" retains spurious "struct" or "union" in the name).
        std::string fullType = m->GetFullTypeName();
        if (fullType != m->GetTrueTypeName()) {
            const std::string& trueName = m->GetTrueTypeName();
            if (fullType.find("::") == std::string::npos && trueName.find("::") != std::string::npos)
                fullType = trueName;
        }

        if ((int)m->GetArrayDim() > 1 || (!m->IsBasic() && m->IsaPointer()))
            fullType.append("*");
        else if ((int)m->GetArrayDim() == 1) {
            std::ostringstream s;
            s << '[' << m->GetMaxIndex(0) << ']' << std::ends;
            fullType.append(s.str());
        }
        return fullType;
    }

    return "<unknown>";
}

intptr_t Cppyy::GetDatamemberOffset(TCppScope_t scope, TCppIndex_t idata)
{
    if (scope == GLOBAL_HANDLE) {
        TGlobal* gbl = g_globalvars[idata];
        if (!gbl->GetAddress() || gbl->GetAddress() == (void*)-1) {
        // CLING WORKAROUND: make sure variable is loaded
            intptr_t addr = (intptr_t)gInterpreter->ProcessLine((std::string("&")+gbl->GetName()+";").c_str());
            if (gbl->GetAddress() && gbl->GetAddress() != (void*)-1)
                return (intptr_t)gbl->GetAddress();        // now loaded!
            return addr;                                   // last resort ...
        }
        return (intptr_t)gbl->GetAddress();
    }

    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass()) {
        TDataMember* m = GetDataMemberByIndex(cr, (int)idata);
    // CLING WORKAROUND: the following causes templates to be instantiated first within the proper
    // scope, making the lookup succeed and preventing spurious duplicate instantiations later. Also,
    // if the variable is not yet loaded, pull it in through gInterpreter.
        if (m->Property() & kIsStatic) {
            if (strchr(cr->GetName(), '<'))
                gInterpreter->ProcessLine(((std::string)cr->GetName()+"::"+m->GetName()+";").c_str());
            if ((intptr_t)m->GetOffsetCint() == (intptr_t)-1)
                return (intptr_t)gInterpreter->ProcessLine((std::string("&")+cr->GetName()+"::"+m->GetName()+";").c_str());
        }
        return (intptr_t)m->GetOffsetCint();    // yes, CINT (GetOffset() is both wrong
                                                // and caches that wrong result!
    }

    return (intptr_t)-1;
}

Cppyy::TCppIndex_t Cppyy::GetDatamemberIndex(TCppScope_t scope, const std::string& name)
{
    if (scope == GLOBAL_HANDLE) {
        TGlobal* gb = (TGlobal*)gROOT->GetListOfGlobals(false /* load */)->FindObject(name.c_str());
        if (!gb) gb = (TGlobal*)gROOT->GetListOfGlobals(true  /* load */)->FindObject(name.c_str());
        if (!gb) {
        // some enums are not loaded as they are not considered part of
        // the global scope, but of the enum scope; get them w/o checking
            TDictionary::DeclId_t did = gInterpreter->GetDataMember(nullptr, name.c_str());
            if (did) {
                DataMemberInfo_t* t = gInterpreter->DataMemberInfo_Factory(did, nullptr);
                ((TListOfDataMembers*)gROOT->GetListOfGlobals())->Get(t, true);
                gb = (TGlobal*)gROOT->GetListOfGlobals(false /* load */)->FindObject(name.c_str());
            }
        }

        if (gb && strcmp(gb->GetFullTypeName(), "(lambda)") == 0) {
        // lambdas use a compiler internal closure type, so we wrap
        // them, then return the wrapper's type
        // TODO: this current leaks the std::function; also, if possible,
        //       should instantiate through TClass rather then ProcessLine
            std::ostringstream s;
            s << "auto __cppyy_internal_wrap_" << name << " = "
                 "new __cling_internal::FT<decltype(" << name << ")>::F"
                 "{" << name << "};";
            gInterpreter->ProcessLine(s.str().c_str());
            TGlobal* wrap = (TGlobal*)gROOT->GetListOfGlobals(true)->FindObject(
                ("__cppyy_internal_wrap_"+name).c_str());
            if (wrap && wrap->GetAddress()) gb = wrap;
        }

        if (gb) {
        // TODO: do we ever need a reverse lookup?
            g_globalvars.push_back(gb);
            return TCppIndex_t(g_globalvars.size() - 1);
        }

    } else {
        TClassRef& cr = type_from_handle(scope);
        if (cr.GetClass()) {
            TDataMember* dm =
                (TDataMember*)cr->GetListOfDataMembers()->FindObject(name.c_str());
            // TODO: turning this into an index is silly ...
            if (dm) return (TCppIndex_t)cr->GetListOfDataMembers()->IndexOf(dm);
            dm = (TDataMember*)cr->GetListOfUsingDataMembers()->FindObject(name.c_str());
            if (dm)
                return (TCppIndex_t)cr->GetListOfDataMembers()->IndexOf(dm)
                    + cr->GetListOfDataMembers()->GetSize();
        }
    }

    return (TCppIndex_t)-1;
}


// data member properties ----------------------------------------------------
bool Cppyy::IsPublicData(TCppScope_t scope, TCppIndex_t idata)
{
    if (scope == GLOBAL_HANDLE)
        return true;
    TClassRef& cr = type_from_handle(scope);
    if (cr->Property() & kIsNamespace)
        return true;
    TDataMember* m = GetDataMemberByIndex(cr, (int)idata);
    return m->Property() & kIsPublic;
}

bool Cppyy::IsProtectedData(TCppScope_t scope, TCppIndex_t idata)
{
    if (scope == GLOBAL_HANDLE)
        return true;
    TClassRef& cr = type_from_handle(scope);
    if (cr->Property() & kIsNamespace)
        return true;
    TDataMember* m = GetDataMemberByIndex(cr, (int)idata);
    return m->Property() & kIsProtected;
}

bool Cppyy::IsStaticData(TCppScope_t scope, TCppIndex_t idata)
{
    if (scope == GLOBAL_HANDLE)
        return true;
    TClassRef& cr = type_from_handle(scope);
    if (cr->Property() & kIsNamespace)
        return true;
    TDataMember* m = GetDataMemberByIndex(cr, (int)idata);
    return m->Property() & kIsStatic;
}

bool Cppyy::IsConstData(TCppScope_t scope, TCppIndex_t idata)
{
    if (scope == GLOBAL_HANDLE) {
        TGlobal* gbl = g_globalvars[idata];
        return gbl->Property() & kIsConstant;
    }
    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass()) {
        TDataMember* m = GetDataMemberByIndex(cr, (int)idata);
        return m->Property() & kIsConstant;
    }
    return false;
}

bool Cppyy::IsEnumData(TCppScope_t scope, TCppIndex_t idata)
{
// TODO: currently, ROOT/meta does not properly distinguish between variables of enum
// type, and values of enums. The latter are supposed to be const. This code relies on
// odd features (bugs?) to figure out the difference, but this should really be fixed
// upstream and/or deserves a new API.

    if (scope == GLOBAL_HANDLE) {
        TGlobal* gbl = g_globalvars[idata];

    // make use of an oddity: enum global variables do not have their kIsStatic bit
    // set, whereas enum global values do
        return (gbl->Property() & kIsEnum) && (gbl->Property() & kIsStatic);
    }

    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass()) {
        TDataMember* m = GetDataMemberByIndex(cr, (int)idata);
        std::string ti = m->GetTypeName();

    // can't check anonymous enums by type name, so just accept them as enums
        if (ti.rfind("(anonymous)") != std::string::npos)
            return m->Property() & kIsEnum;

    // since there seems to be no distinction between data of enum type and enum values,
    // check the list of constants for the type to see if there's a match
        if (ti.rfind(cr->GetName(), 0) != std::string::npos) {
            std::string::size_type s = strlen(cr->GetName())+2;
            if (s < ti.size()) {
                TEnum* ee = ((TListOfEnums*)cr->GetListOfEnums())->GetObject(ti.substr(s, std::string::npos).c_str());
                if (ee) return ee->GetConstant(m->GetName());
            }
        }
    }

// this default return only means that the data will be writable, not that it will
// be unreadable or otherwise misrepresented
    return false;
}

int Cppyy::GetDimensionSize(TCppScope_t scope, TCppIndex_t idata, int dimension)
{
    if (scope == GLOBAL_HANDLE) {
        TGlobal* gbl = g_globalvars[idata];
        return gbl->GetMaxIndex(dimension);
    }
    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass()) {
        TDataMember* m = GetDataMemberByIndex(cr, (int)idata);
        return m->GetMaxIndex(dimension);
    }
    return -1;
}


// enum properties -----------------------------------------------------------
Cppyy::TCppEnum_t Cppyy::GetEnum(TCppScope_t scope, const std::string& enum_name)
{
    if (scope == GLOBAL_HANDLE)
        return (TCppEnum_t)gROOT->GetListOfEnums(kTRUE)->FindObject(enum_name.c_str());

    TClassRef& cr = type_from_handle(scope);
    if (cr.GetClass())
        return (TCppEnum_t)cr->GetListOfEnums(kTRUE)->FindObject(enum_name.c_str());

    return (TCppEnum_t)0;
}

Cppyy::TCppIndex_t Cppyy::GetNumEnumData(TCppEnum_t etype)
{
    return (TCppIndex_t)((TEnum*)etype)->GetConstants()->GetSize();
}

std::string Cppyy::GetEnumDataName(TCppEnum_t etype, TCppIndex_t idata)
{
    return ((TEnumConstant*)((TEnum*)etype)->GetConstants()->At(idata))->GetName();
}

long long Cppyy::GetEnumDataValue(TCppEnum_t etype, TCppIndex_t idata)
{
     TEnumConstant* ecst = (TEnumConstant*)((TEnum*)etype)->GetConstants()->At(idata);
     return (long long)ecst->GetValue();
}


//- C-linkage wrappers -------------------------------------------------------

extern "C" {
/* direct interpreter access ---------------------------------------------- */
int cppyy_compile(const char* code) {
    return Cppyy::Compile(code);
}


/* name to opaque C++ scope representation -------------------------------- */
char* cppyy_resolve_name(const char* cppitem_name) {
    return cppstring_to_cstring(Cppyy::ResolveName(cppitem_name));
}

char* cppyy_resolve_enum(const char* enum_type) {
    return cppstring_to_cstring(Cppyy::ResolveEnum(enum_type));
}

cppyy_scope_t cppyy_get_scope(const char* scope_name) {
    return cppyy_scope_t(Cppyy::GetScope(scope_name));
}

cppyy_type_t cppyy_actual_class(cppyy_type_t klass, cppyy_object_t obj) {
    return cppyy_type_t(Cppyy::GetActualClass(klass, (void*)obj));
}

size_t cppyy_size_of_klass(cppyy_type_t klass) {
    return Cppyy::SizeOf(klass);
}

size_t cppyy_size_of_type(const char* type_name) {
    return Cppyy::SizeOf(type_name);
}


/* memory management ------------------------------------------------------ */
cppyy_object_t cppyy_allocate(cppyy_type_t type) {
    return cppyy_object_t(Cppyy::Allocate(type));
}

void cppyy_deallocate(cppyy_type_t type, cppyy_object_t self) {
    Cppyy::Deallocate(type, (void*)self);
}

cppyy_object_t cppyy_construct(cppyy_type_t type) {
    return (cppyy_object_t)Cppyy::Construct(type);
}

void cppyy_destruct(cppyy_type_t type, cppyy_object_t self) {
    Cppyy::Destruct(type, (void*)self);
}


/* method/function dispatching -------------------------------------------- */
/* Exception types:
    1: default (unknown exception)
    2: standard exception
*/
#define CPPYY_HANDLE_EXCEPTION                                               \
    catch (std::exception& e) {                                              \
        cppyy_exctype_t* etype = (cppyy_exctype_t*)((Parameter*)args+nargs); \
        *etype = (cppyy_exctype_t)2;                                         \
        *((char**)(etype+1)) = cppstring_to_cstring(e.what());               \
    }                                                                        \
    catch (...) {                                                            \
        cppyy_exctype_t* etype = (cppyy_exctype_t*)((Parameter*)args+nargs); \
        *etype = (cppyy_exctype_t)1;                                         \
        *((char**)(etype+1)) =                                               \
            cppstring_to_cstring("unhandled, unknown C++ exception");        \
    }

void cppyy_call_v(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    try {
        Cppyy::CallV(method, (void*)self, nargs, args);
    } CPPYY_HANDLE_EXCEPTION
}

unsigned char cppyy_call_b(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    try {
        return (unsigned char)Cppyy::CallB(method, (void*)self, nargs, args);
    } CPPYY_HANDLE_EXCEPTION
    return (unsigned char)-1;
}

char cppyy_call_c(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    try {
        return (char)Cppyy::CallC(method, (void*)self, nargs, args);
    } CPPYY_HANDLE_EXCEPTION
    return (char)-1;
}

short cppyy_call_h(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    try {
        return (short)Cppyy::CallH(method, (void*)self, nargs, args);
    } CPPYY_HANDLE_EXCEPTION
    return (short)-1;
}

int cppyy_call_i(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    try {
        return (int)Cppyy::CallI(method, (void*)self, nargs, args);
    } CPPYY_HANDLE_EXCEPTION
    return (int)-1;
}

long cppyy_call_l(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    try {
        return (long)Cppyy::CallL(method, (void*)self, nargs, args);
    } CPPYY_HANDLE_EXCEPTION
    return (long)-1;
}

long long cppyy_call_ll(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    try {
        return (long long)Cppyy::CallLL(method, (void*)self, nargs, args);
    } CPPYY_HANDLE_EXCEPTION
    return (long long)-1;
}

float cppyy_call_f(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    try {
        return (float)Cppyy::CallF(method, (void*)self, nargs, args);
    } CPPYY_HANDLE_EXCEPTION
    return (float)-1;
}

double cppyy_call_d(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    try {
        return (double)Cppyy::CallD(method, (void*)self, nargs, args);
    } CPPYY_HANDLE_EXCEPTION
    return (double)-1;
}

long double cppyy_call_ld(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    try {
        return (long double)Cppyy::CallLD(method, (void*)self, nargs, args);
    } CPPYY_HANDLE_EXCEPTION
    return (long double)-1;
}

double cppyy_call_nld(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    return (double)cppyy_call_ld(method, self, nargs, args);
}

void* cppyy_call_r(cppyy_method_t method, cppyy_object_t self, int nargs, void* args) {
    try {
        return (void*)Cppyy::CallR(method, (void*)self, nargs, args);
    } CPPYY_HANDLE_EXCEPTION
    return (void*)nullptr;
}

char* cppyy_call_s(
        cppyy_method_t method, cppyy_object_t self, int nargs, void* args, size_t* lsz) {
    try {
        return Cppyy::CallS(method, (void*)self, nargs, args, lsz);
    } CPPYY_HANDLE_EXCEPTION
    return (char*)nullptr;
}

cppyy_object_t cppyy_constructor(
        cppyy_method_t method, cppyy_type_t klass, int nargs, void* args) {
    try {
        return cppyy_object_t(Cppyy::CallConstructor(method, klass, nargs, args));
    } CPPYY_HANDLE_EXCEPTION
    return (cppyy_object_t)0;
}

void cppyy_destructor(cppyy_type_t klass, cppyy_object_t self) {
    Cppyy::CallDestructor(klass, self);
}

cppyy_object_t cppyy_call_o(cppyy_method_t method, cppyy_object_t self,
        int nargs, void* args, cppyy_type_t result_type) {
    try {
        return cppyy_object_t(Cppyy::CallO(method, (void*)self, nargs, args, result_type));
    } CPPYY_HANDLE_EXCEPTION
    return (cppyy_object_t)0;
}

cppyy_funcaddr_t cppyy_function_address(cppyy_method_t method) {
    return cppyy_funcaddr_t(Cppyy::GetFunctionAddress(method, true));
}


/* handling of function argument buffer ----------------------------------- */
void* cppyy_allocate_function_args(int nargs) {
// for calls through C interface, require extra space for reporting exceptions
    return malloc(nargs*sizeof(Parameter)+sizeof(cppyy_exctype_t)+sizeof(char**));
}

void cppyy_deallocate_function_args(void* args) {
    free(args);
}

size_t cppyy_function_arg_sizeof() {
    return (size_t)Cppyy::GetFunctionArgSizeof();
}

size_t cppyy_function_arg_typeoffset() {
    return (size_t)Cppyy::GetFunctionArgTypeoffset();
}


/* scope reflection information ------------------------------------------- */
int cppyy_is_namespace(cppyy_scope_t scope) {
    return (int)Cppyy::IsNamespace(scope);
}

int cppyy_is_template(const char* template_name) {
    return (int)Cppyy::IsTemplate(template_name);
}

int cppyy_is_abstract(cppyy_type_t type) {
    return (int)Cppyy::IsAbstract(type);
}

int cppyy_is_enum(const char* type_name) {
    return (int)Cppyy::IsEnum(type_name);
}

const char** cppyy_get_all_cpp_names(cppyy_scope_t scope, size_t* count) {
    std::set<std::string> cppnames;
    Cppyy::GetAllCppNames(scope, cppnames);
    const char** c_cppnames = (const char**)malloc(cppnames.size()*sizeof(const char*));
    int i = 0;
    for (const auto& name : cppnames) {
        c_cppnames[i] = cppstring_to_cstring(name);
        ++i;
    }
    *count = cppnames.size();
    return c_cppnames;
}


/* namespace reflection information --------------------------------------- */
cppyy_scope_t* cppyy_get_using_namespaces(cppyy_scope_t scope) {
    const std::vector<Cppyy::TCppScope_t>& uv = Cppyy::GetUsingNamespaces((Cppyy::TCppScope_t)scope);

    if (uv.empty())
        return (cppyy_index_t*)nullptr;

    cppyy_scope_t* llresult = (cppyy_scope_t*)malloc(sizeof(cppyy_scope_t)*(uv.size()+1));
    for (int i = 0; i < (int)uv.size(); ++i) llresult[i] = uv[i];
    llresult[uv.size()] = (cppyy_scope_t)0;
    return llresult;
}


/* class reflection information ------------------------------------------- */
char* cppyy_final_name(cppyy_type_t type) {
    return cppstring_to_cstring(Cppyy::GetFinalName(type));
}

char* cppyy_scoped_final_name(cppyy_type_t type) {
    return cppstring_to_cstring(Cppyy::GetScopedFinalName(type));
}

int cppyy_has_virtual_destructor(cppyy_type_t type) {
    return (int)Cppyy::HasVirtualDestructor(type);
}

int cppyy_has_complex_hierarchy(cppyy_type_t type) {
    return (int)Cppyy::HasComplexHierarchy(type);
}

int cppyy_num_bases(cppyy_type_t type) {
    return (int)Cppyy::GetNumBases(type);
}

char* cppyy_base_name(cppyy_type_t type, int base_index) {
    return cppstring_to_cstring(Cppyy::GetBaseName (type, base_index));
}

int cppyy_is_subtype(cppyy_type_t derived, cppyy_type_t base) {
    return (int)Cppyy::IsSubtype(derived, base);
}

int cppyy_is_smartptr(cppyy_type_t type) {
    return (int)Cppyy::IsSmartPtr(type);
}

int cppyy_smartptr_info(const char* name, cppyy_type_t* raw, cppyy_method_t* deref) {
    return (int)Cppyy::GetSmartPtrInfo(name, raw, deref);
}

void cppyy_add_smartptr_type(const char* type_name) {
    Cppyy::AddSmartPtrType(type_name);
}


/* calculate offsets between declared and actual type, up-cast: direction > 0; down-cast: direction < 0 */
ptrdiff_t cppyy_base_offset(cppyy_type_t derived, cppyy_type_t base, cppyy_object_t address, int direction) {
    return (ptrdiff_t)Cppyy::GetBaseOffset(derived, base, (void*)address, direction, 0);
}


/* method/function reflection information --------------------------------- */
int cppyy_num_methods(cppyy_scope_t scope) {
    return (int)Cppyy::GetNumMethods(scope);
}

cppyy_index_t* cppyy_method_indices_from_name(cppyy_scope_t scope, const char* name)
{
    std::vector<cppyy_index_t> result = Cppyy::GetMethodIndicesFromName(scope, name);

    if (result.empty())
        return (cppyy_index_t*)nullptr;

    cppyy_index_t* llresult = (cppyy_index_t*)malloc(sizeof(cppyy_index_t)*(result.size()+1));
    for (int i = 0; i < (int)result.size(); ++i) llresult[i] = result[i];
    llresult[result.size()] = -1;
    return llresult;
}

cppyy_method_t cppyy_get_method(cppyy_scope_t scope, cppyy_index_t idx) {
    return cppyy_method_t(Cppyy::GetMethod(scope, idx));
}

char* cppyy_method_name(cppyy_method_t method) {
    return cppstring_to_cstring(Cppyy::GetMethodName((Cppyy::TCppMethod_t)method));
}

char* cppyy_method_full_name(cppyy_method_t method) {
    return cppstring_to_cstring(Cppyy::GetMethodFullName((Cppyy::TCppMethod_t)method));
}

char* cppyy_method_mangled_name(cppyy_method_t method) {
    return cppstring_to_cstring(Cppyy::GetMethodMangledName((Cppyy::TCppMethod_t)method));
}

char* cppyy_method_result_type(cppyy_method_t method) {
    return cppstring_to_cstring(Cppyy::GetMethodResultType((Cppyy::TCppMethod_t)method));
}

int cppyy_method_num_args(cppyy_method_t method) {
    return (int)Cppyy::GetMethodNumArgs((Cppyy::TCppMethod_t)method);
}

int cppyy_method_req_args(cppyy_method_t method) {
    return (int)Cppyy::GetMethodReqArgs((Cppyy::TCppMethod_t)method);
}

char* cppyy_method_arg_name(cppyy_method_t method, int arg_index) {
    return cppstring_to_cstring(Cppyy::GetMethodArgName((Cppyy::TCppMethod_t)method, (Cppyy::TCppIndex_t)arg_index));
}

char* cppyy_method_arg_type(cppyy_method_t method, int arg_index) {
    return cppstring_to_cstring(Cppyy::GetMethodArgType((Cppyy::TCppMethod_t)method, (Cppyy::TCppIndex_t)arg_index));
}

char* cppyy_method_arg_default(cppyy_method_t method, int arg_index) {
    return cppstring_to_cstring(Cppyy::GetMethodArgDefault((Cppyy::TCppMethod_t)method, (Cppyy::TCppIndex_t)arg_index));
}

char* cppyy_method_signature(cppyy_method_t method, int show_formalargs) {
    return cppstring_to_cstring(Cppyy::GetMethodSignature((Cppyy::TCppMethod_t)method, (bool)show_formalargs));
}

char* cppyy_method_signature_max(cppyy_method_t method, int show_formalargs, int maxargs) {
    return cppstring_to_cstring(Cppyy::GetMethodSignature((Cppyy::TCppMethod_t)method, (bool)show_formalargs, (Cppyy::TCppIndex_t)maxargs));
}

char* cppyy_method_prototype(cppyy_scope_t scope, cppyy_method_t method, int show_formalargs) {
    return cppstring_to_cstring(Cppyy::GetMethodPrototype(
        (Cppyy::TCppScope_t)scope, (Cppyy::TCppMethod_t)method, (bool)show_formalargs));
}

int cppyy_is_const_method(cppyy_method_t method) {
    return (int)Cppyy::IsConstMethod(method);
}

int cppyy_get_num_templated_methods(cppyy_scope_t scope) {
    return (int)Cppyy::GetNumTemplatedMethods(scope);
}

char* cppyy_get_templated_method_name(cppyy_scope_t scope, cppyy_index_t imeth) {
    return cppstring_to_cstring(Cppyy::GetTemplatedMethodName(scope, imeth));
}

int cppyy_is_templated_constructor(cppyy_scope_t scope, cppyy_index_t imeth) {
    return Cppyy::IsTemplatedConstructor((Cppyy::TCppScope_t)scope, (Cppyy::TCppIndex_t)imeth);
}

int cppyy_exists_method_template(cppyy_scope_t scope, const char* name) {
    return (int)Cppyy::ExistsMethodTemplate(scope, name);
}

int cppyy_method_is_template(cppyy_scope_t scope, cppyy_index_t idx) {
    return (int)Cppyy::IsMethodTemplate(scope, idx);
}

cppyy_method_t cppyy_get_method_template(cppyy_scope_t scope, const char* name, const char* proto) {
    return cppyy_method_t(Cppyy::GetMethodTemplate(scope, name, proto));
}

cppyy_index_t cppyy_get_global_operator(cppyy_scope_t scope, cppyy_scope_t lc, cppyy_scope_t rc, const char* op) {
    return cppyy_index_t(Cppyy::GetGlobalOperator(scope, Cppyy::GetScopedFinalName(lc), Cppyy::GetScopedFinalName(rc), op));
}


/* method properties ------------------------------------------------------ */
int cppyy_is_publicmethod(cppyy_method_t method) {
    return (int)Cppyy::IsPublicMethod((Cppyy::TCppMethod_t)method);
}

int cppyy_is_protectedmethod(cppyy_method_t method) {
    return (int)Cppyy::IsProtectedMethod((Cppyy::TCppMethod_t)method);
}

int cppyy_is_constructor(cppyy_method_t method) {
    return (int)Cppyy::IsConstructor((Cppyy::TCppMethod_t)method);
}

int cppyy_is_destructor(cppyy_method_t method) {
    return (int)Cppyy::IsDestructor((Cppyy::TCppMethod_t)method);
}

int cppyy_is_staticmethod(cppyy_method_t method) {
    return (int)Cppyy::IsStaticMethod((Cppyy::TCppMethod_t)method);
}


/* data member reflection information ------------------------------------- */
int cppyy_num_datamembers(cppyy_scope_t scope) {
    return (int)Cppyy::GetNumDatamembers(scope);
}

char* cppyy_datamember_name(cppyy_scope_t scope, int datamember_index) {
    return cppstring_to_cstring(Cppyy::GetDatamemberName(scope, datamember_index));
}

char* cppyy_datamember_type(cppyy_scope_t scope, int datamember_index) {
    return cppstring_to_cstring(Cppyy::GetDatamemberType(scope, datamember_index));
}

intptr_t cppyy_datamember_offset(cppyy_scope_t scope, int datamember_index) {
    return intptr_t(Cppyy::GetDatamemberOffset(scope, datamember_index));
}

int cppyy_datamember_index(cppyy_scope_t scope, const char* name) {
    return (int)Cppyy::GetDatamemberIndex(scope, name);
}



/* data member properties ------------------------------------------------- */
int cppyy_is_publicdata(cppyy_type_t type, cppyy_index_t datamember_index) {
    return (int)Cppyy::IsPublicData(type, datamember_index);
}

int cppyy_is_protecteddata(cppyy_type_t type, cppyy_index_t datamember_index) {
    return (int)Cppyy::IsProtectedData(type, datamember_index);
}

int cppyy_is_staticdata(cppyy_type_t type, cppyy_index_t datamember_index) {
    return (int)Cppyy::IsStaticData(type, datamember_index);
}

int cppyy_is_const_data(cppyy_scope_t scope, cppyy_index_t idata) {
    return (int)Cppyy::IsConstData(scope, idata);
}

int cppyy_is_enum_data(cppyy_scope_t scope, cppyy_index_t idata) {
    return (int)Cppyy::IsEnumData(scope, idata);
}

int cppyy_get_dimension_size(cppyy_scope_t scope, cppyy_index_t idata, int dimension) {
    return Cppyy::GetDimensionSize(scope, idata, dimension);
}


/* misc helpers ----------------------------------------------------------- */
RPY_EXTERN
void* cppyy_load_dictionary(const char* lib_name) {
    int result = gSystem->Load(lib_name);
    return (void*)(result == 0 /* success */ || result == 1 /* already loaded */);
}

#if defined(_MSC_VER)
long long cppyy_strtoll(const char* str) {
    return _strtoi64(str, NULL, 0);
}

extern "C" {
unsigned long long cppyy_strtoull(const char* str) {
    return _strtoui64(str, NULL, 0);
}
}
#else
long long cppyy_strtoll(const char* str) {
    return strtoll(str, NULL, 0);
}

extern "C" {
unsigned long long cppyy_strtoull(const char* str) {
    return strtoull(str, NULL, 0);
}
}
#endif

void cppyy_free(void* ptr) {
    free(ptr);
}

cppyy_object_t cppyy_charp2stdstring(const char* str, size_t sz) {
    return (cppyy_object_t)new std::string(str, sz);
}

const char* cppyy_stdstring2charp(cppyy_object_t ptr, size_t* lsz) {
    *lsz = ((std::string*)ptr)->size();
    return ((std::string*)ptr)->data();
}

cppyy_object_t cppyy_stdstring2stdstring(cppyy_object_t ptr) {
    return (cppyy_object_t)new std::string(*(std::string*)ptr);
}

double cppyy_longdouble2double(void* p) {
    return (double)*(long double*)p;
}

void cppyy_double2longdouble(double d, void* p) {
    *(long double*)p = d;
}

int cppyy_vectorbool_getitem(cppyy_object_t ptr, int idx) {
    return (int)(*(std::vector<bool>*)ptr)[idx];
}

void cppyy_vectorbool_setitem(cppyy_object_t ptr, int idx, int value) {
    (*(std::vector<bool>*)ptr)[idx] = (bool)value;
}

} // end C-linkage wrappers
