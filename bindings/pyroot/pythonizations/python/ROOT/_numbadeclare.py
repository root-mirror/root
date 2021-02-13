from cppyy import gbl as gbl_namespace


def _NumbaDeclareDecorator(input_types, return_type, name=None):
    '''
    Decorator for making Python callables accessible in C++ by just-in-time compilation
    with numba and cling

    The decorator takes the given Python callable and just-in-time compiles (jits)
    wrapper functions with the given C++ types for input and return types. Eventually,
    the Python callable is accessible in the Numba namespace in C++.

    The implementation first jits with numba the Python callable. We support fundamental types and
    ROOT::VecOps::RVecs thereof. Note that you can get the jitted Python callable by the attribute
    numba_func. The C++ types are converted to the respective numba types and RVecs are accessible
    in Python by numpy arrays. After jitting the actual Python callable, we jit another Python wrapper,
    which converts the Python signature to a C-friendly signature. The wrapper code is accessible by
    the attribute __py_wrapper__. Next, the Python wrapper is given to cling to jit a C++ wrapper function,
    making the original Python callable accessible in C++. The wrapper code in C++ is accessible by
    the attribute __cpp_wrapper__.

    Note that the callable is fully compiled without side-effects. The numba jitting uses the nopython
    option which does not allow interaction with the Python interpreter. This means that you can use
    the resulting function also safely in multi-threaded environments.
    '''
    # Make required imports
    try:
        import numba as nb
    except:
        raise Exception('Failed to import numba')
    try:
        import cffi
    except:
        raise Exception('Failed to import cffi')
    import re, sys

    # Normalize input types by stripping ROOT and VecOps namespaces from input types
    def normalize_typename(t):
        '''
        Remove ROOT:: and VecOps:: namespaces
        '''
        return t.replace('ROOT::', '').replace('VecOps::', '')

    input_types = [normalize_typename(t) for t in input_types]
    return_type = normalize_typename(return_type)

    # Helper functions to determine types
    def get_inner_type(t):
        '''
        Get inner typename of a templated C++ typename
        '''
        try:
            g = re.match('(.*)<(.*)>', t).groups(0)
            return g[1]
        except:
            raise Exception('Failed to extract template argument of type {}'.format(t))

    def get_numba_type(t):
        '''
        Get numba type object from a C++ fundamental typename

        These are the types we use to jit the Python callable.
        '''
        typemap = {
                'float': nb.float32,
                'double': nb.float64,
                'int': nb.int32,
                'unsigned int': nb.uint32,
                'long': nb.int64,
                'unsigned long': nb.uint64,
                'bool': nb.boolean
                }
        if t in typemap:
            return typemap[t]
        raise Exception(
                'Type {} is not supported for jitting with numba. Valid fundamental types and RVecs thereof are {}'.format(
                    t, list(typemap.keys())))

    def get_c_signature(input_types, return_type):
        '''
        Get C friendly signature as numba type objects from C++ typenames

        We need the types to jit a Python wrapper, which can be accessed as a function pointer in C++.
        '''
        c_input_types = []
        for t in input_types:
            if 'RVec' in t:
                c_input_types += [nb.types.CPointer(get_numba_type(get_inner_type(t))), nb.int32]
            else:
                c_input_types.append(get_numba_type(t))
        if 'RVec' in return_type:
            # We return an RVec through pointers as part of the input arguments. Note that the
            # pointer type in numba is always an int64 and is later on cast in C++ to the correct type.
            # In addition, we provide the size of the data type of the array for preallocating memory of
            # the returned array.
            # See the Python wrapper for further information why we are using these types.
            c_return_type = nb.void
            c_input_types += [
                    nb.types.CPointer(nb.int64), # Pointer to the data (the first element of the array)
                    nb.types.CPointer(nb.int64)] # Size of the array in elements
        else:
            c_return_type = get_numba_type(return_type)
        return c_return_type, c_input_types

    def get_numba_signature(input_types, return_type):
        '''
        Get numba signature as numba type objects from C++ typenames
        '''
        nb_input_types = []
        for t in input_types:
            if 'RVec' in t:
                nb_input_types.append(get_numba_type(get_inner_type(t))[:])
            else:
                nb_input_types.append(get_numba_type(t))
        if 'RVec' in return_type:
            nb_return_type = get_numba_type(get_inner_type(return_type))[:]
        else:
            nb_return_type = get_numba_type(return_type)
        return nb_return_type, nb_input_types

    def inner(func, input_types=input_types, return_type=return_type, name=name):
        '''
        Inner decorator without arguments, see outer decorator for documentation
        '''

        # Jit the given Python callable with numba
        nb_return_type, nb_input_types = get_numba_signature(input_types, return_type)
        try:
            nbjit = nb.jit(nb_return_type(*nb_input_types), nopython=True, inline='always')(func)
        except:
            raise Exception('Failed to jit Python callable {} with numba.jit'.format(func))
        func.numba_func = nbjit

        # Create Python wrapper with C++ friendly signature

        # Define signature
        pywrapper_signature = [
                'ptr_{0}, size_{0}'.format(i) if 'RVec' in t else 'x_{}'.format(i) \
                        for i, t in enumerate(input_types)]
        if 'RVec' in return_type:
            # If we return an RVec, we return via pointer the pointer of the allocated data,
            # the size in elements. In addition, we provide the size of the datatype in bytes.
            pywrapper_signature += ['ptrptr_r, ptrsize_r']

        # Define arguments for jit function
        pywrapper_args_def = [
                'x_{0} = nb.carray(ptr_{0}, (size_{0},))'.format(i) if 'RVec' in t else 'x_{}'.format(i) \
                        for i, t in enumerate(input_types)]
        pywrapper_args = ['x_{}'.format(i) for i in range(len(input_types))]

        # Define return operation
        if 'RVec' in return_type:
            innert = get_inner_type(return_type)
            dtypesize = 1 if innert == 'bool' else int(get_numba_type(innert).bitwidth / 8)
            pywrapper_return = '\n    '.join([
                '# Because we cannot manipulate the memory management of the numpy array we copy the data',
                'ptr = malloc(r.size * {})'.format(dtypesize),
                'cp = nb.carray(ptr, r.size, dtype_r)',
                'cp[:] = r[:]',
                '# Return size of the array and the pointer to the copied data',
                'ptrsize_r[0] = r.size',
                'ptrptr_r[0] = cp.ctypes.data'
                ])
        else:
            pywrapper_return = 'return r'

        # Build wrapper code
        pywrappercode = '''\
def pywrapper({SIGNATURE}):
    """
    Wrapper function for the jitted Python callable with special treatment of arrays
    """
    # If an RVec is given, define numba carray wrapper for the input types
    {ARGS_DEF}
    # Call the jitted Python function
    r = nbjit({ARGS})
    # Return the result
    {RETURN}
        '''.format(
                SIGNATURE=', '.join(pywrapper_signature),
                ARGS_DEF='\n    '.join(pywrapper_args_def),
                ARGS=', '.join(pywrapper_args),
                RETURN=pywrapper_return
                )

        glob = dict(globals()) # Make a shallow copy of the dictionary so we don't pollute the global scope
        glob['nb'] = nb
        glob['nbjit'] = nbjit

        ffi = cffi.FFI()
        ffi.cdef('void* malloc(long size);')
        C = ffi.dlopen(None)
        glob['malloc'] = C.malloc

        if 'RVec' in return_type:
            glob['dtype_r'] = get_numba_type(get_inner_type(return_type))

        if sys.version_info[0] >= 3:
            exec(pywrappercode, glob, locals()) in {}
        else:
            exec(pywrappercode) in glob, locals()

        if not 'pywrapper' in locals():
            raise Exception('Failed to create Python wrapper function:\n{}'.format(pywrappercode))

        # Jit the Python wrapper code
        c_return_type, c_input_types = get_c_signature(input_types, return_type)
        try:
            nbcfunc = nb.cfunc(c_return_type(*c_input_types), nopython=True)(locals()['pywrapper'])
        except:
            raise Exception('Failed to jit Python wrapper with numba.cfunc')
        func.__py_wrapper__ = pywrappercode
        func.__numba_cfunc__ = nbcfunc

        # Get address of jitted wrapper function
        address = nbcfunc.address

        # Infer name of the C++ wrapper function
        if not name:
            name = func.__name__

        # Build C++ wrapper for jitting with cling

        # Define input signature
        input_types_ref = ['ROOT::{}&'.format(t) if 'RVec' in t else t for t in input_types]
        input_signature = ', '.join('{} x_{}'.format(t, i) for i, t in enumerate(input_types_ref))

        # Define function pointer types
        func_ptr_input_types = []
        for t in input_types:
            if 'RVec' in t:
                innert = get_inner_type(t)
                if innert == 'bool':
                    # Special treatment for bool: In numpy, bools have 1 byte
                    func_ptr_input_types += ['char*, int']
                else:
                    func_ptr_input_types += ['{}*, int'.format(innert)]
            else:
                func_ptr_input_types += [t]
        if 'RVec' in return_type:
            # See C++ wrapper code for the reason using these types
            innert = get_inner_type(return_type)
            func_ptr_input_types += ['{}**, long*'.format('char' if innert == 'bool' else innert)]
        func_ptr_type = '{RETURN_TYPE}(*)({INPUT_TYPES})'.format(
                RETURN_TYPE='void*' if 'RVec' in return_type else return_type,
                INPUT_TYPES=', '.join(func_ptr_input_types)
                )

        # Define function call
        vecbool_conversion = []
        func_args = []
        for i, t in enumerate(input_types):
            if 'RVec' in t:
                func_args += ['x_{0}.data(), x_{0}.size()'.format(i)]
                if get_inner_type(t) == 'bool':
                    # Copy the RVec<bool> to a RVec<char> to match the numpy memory layout
                    func_args[-1] = func_args[-1].replace('x_', 'xb_')
                    vecbool_conversion += ['ROOT::RVec<char> xb_{0} = x_{0};'.format(i)]
            else:
                func_args += ['x_{}'.format(i)]
        if 'RVec' in return_type:
            # See C++ wrapper code for the reason using these arguments
            func_args += ['&ptr, &size']

        # Define return operation
        if 'RVec' in return_type:
            innert = get_inner_type(return_type)
            if innert == 'bool': innert = 'char'
            return_op = '\n    '.join([
                '// Because an RVec cannot take the ownership of external data, we have to copy the returned array',
                'long size; // Size of the returned array',
                '{}* ptr; // Pointer to the data of the returned array'.format(innert),
                'funcptr({});'.format(', '.join(func_args)),
                # TODO: Remove this copy as soon as RVec can adopt the ownership
                'ROOT::RVec<{}> x_r(ptr, ptr + size);'.format(innert),
                'free(ptr);',
                # If we return a RVec<bool>, we rely here on the automatic conversion of RVec<char> to RVec<bool>
                'return x_r;'])
        else:
            return_op = 'return funcptr({});'.format(', '.join(func_args))

        # Build wrapper code
        cppwrappercode = """\
namespace Numba {{
/*
 * C++ wrapper function around the jitted Python wrapping which calls the jitted Python callable
 */
{RETURN_TYPE} {FUNC_NAME}({INPUT_SIGNATURE}) {{
    // Create a function pointer from the jitted Python wrapper
    const auto funcptr = reinterpret_cast<{FUNC_PTR_TYPE}>({FUNC_PTR});
    // Perform conversion of RVec<bool>
    {VECBOOL_CONVERSION}
    // Return the result
    {RETURN_OP}
}}
}}""".format(
                RETURN_TYPE='ROOT::' + return_type if 'RVec' in return_type else return_type,
                FUNC_NAME=name,
                INPUT_SIGNATURE=input_signature,
                FUNC_PTR=address,
                FUNC_PTR_TYPE=func_ptr_type,
                VECBOOL_CONVERSION='\n    '.join(vecbool_conversion),
                RETURN_OP=return_op)

        # Jit wrapper C++ code
        err = gbl_namespace.gInterpreter.Declare(cppwrappercode)
        if not err:
            raise Exception('Failed to jit C++ wrapper code with cling:\n{}'.format(cppwrappercode))
        func.__cpp_wrapper__ = cppwrappercode

        return func

    return inner
