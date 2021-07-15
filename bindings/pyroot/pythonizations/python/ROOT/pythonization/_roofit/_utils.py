# Authors:
# * Hinnerk C. Schmidt 02/2021
# * Jonas Rembser 06/2021
# * Harshal Shende 06/2021

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################


def _getter(k, v):
    # helper function to get CmdArg attribute from `RooFit`
    # Parameters:
    # k: key of the kwarg
    # v: value of the kwarg

    # We have to use ROOT here and not cppy.gbl, because the RooFit namespace is pythonized itself.
    import ROOT
    import libcppyy

    func = getattr(ROOT.RooFit, k)

    if isinstance(func, libcppyy.CPPOverload):
        # Pythonization for functions that don't pass any RooCmdArgs like ShiftToZero() and MoveToBack(). For Eg,
        # Default bindings: pdf.plotOn(frame, ROOT.RooFit.MoveToBack())
        # With pythonizations: pdf.plotOn(frame, MoveToBack=True)

        if "()" in func.func_doc:
            if not isinstance(v, bool):
                raise TypeError("The keyword argument " + k + " can only take bool values.")
            return func() if v else ROOT.RooCmdArg.none()

    try:
        # If the keyword argument value is a tuple, list, or dict, first
        # try to unpack it as parameters to the RooCmdArg-generating
        # function. If this doesn't succeed, the tuple, list, or dict,
        # will be passed directly to the function as it's only argument.
        if isinstance(v, (tuple, list)):
            return func(*v)
        elif isinstance(v, (dict,)):
            return func(**v)
    except:
        pass

    return func(v)


def _kwargs_to_roocmdargs(*args, **kwargs):
    """Helper function to check kwargs and pythonize the arguments using _getter"""
    if kwargs:
        args = args + tuple((_getter(k, v) for k, v in kwargs.items()))
    return args, {}


def _string_to_root_attribute(value, lookup_map):
    """Helper function to pythonize arguments based on the matplotlib color/style conventions."""
    # lookup_map for color and style defined in _rooglobalfunc.py have matplotlib conventions specified which enables
    # the use of string values like "r" or ":" instead of using enums like ROOT.kRed or ROOT.kDotted for colors or styles. For Eg.
    # Default bindings:  pdf.plotOn(frame, LineColor=ROOT.kOrange)
    # With pythonizations: pdf.plotOn(frame, LineColor="kOrange")

    import ROOT

    if isinstance(value, str):
        if value in lookup_map:
            return getattr(ROOT, lookup_map[value])
        else:
            try:
                return getattr(ROOT, value)
            except:
                raise ValueError(
                    "Unsupported value passed. The value either has to be the name of an attribute of the ROOT module, or match with one of the following values that get translated to ROOT attributes: {}".format(
                        _lookup_map
                    )
                )
    else:
        return value


def _dict_to_std_map(arg_dict):
    """Helper function to convert python dict to std::map."""
    import ROOT

    allowed_key_types = ["std::string", "RooCategory*"]
    allowed_value_types = ["RooDataSet*", "RooDataHist*", "TH1*", "std::string", "int"]

    def all_of_class(d, type, check_key):
        return all([isinstance(key if check_key else value, type) for key, value in d.items()])

    def get_python_class(cpp_type_name):

        cpp_type_name = cpp_type_name.replace("*", "")

        if cpp_type_name in ["std::string", "string"]:
            return str
        if cpp_type_name == "int":
            return int

        # otherwise try to get class from the ROOT namespace
        return getattr(ROOT, cpp_type_name)

    def get_template_args(import_dict):

        key_type = None
        value_type = None

        for typename in allowed_key_types:
            if all_of_class(import_dict, get_python_class(typename), True):
                key_type = typename

        for typename in allowed_value_types:
            if all_of_class(import_dict, get_python_class(typename), False):
                value_type = typename

        def get_python_typenames(typenames):
            return [get_python_class(t).__name__ for t in typenames]

        def prettyprint_str_list(l):
            if len(l) == 1:
                return l[0]
            if len(l) == 2:
                return l[0] + " or " + l[1]
            return ", ".join(l[:-1]) + ", or " + l[-1]

        if key_type is None:
            raise TypeError(
                "All dictionary keys must be of the same type, which can be either "
                + prettyprint_str_list(get_python_typenames(allowed_key_types))
                + "."
            )

        if value_type is None:
            raise TypeError(
                "All dictionary values must be of the same type, which can be either "
                + prettyprint_str_list(get_python_typenames(allowed_value_types))
                + "."
            )

        return key_type + "," + value_type

    # The std::map created by this function usually contains non-owning pointers as values.
    # This is not considered by Pythons reference counter. To ensure that the pointed-to objects
    # live at least as long as the std::map, a python list containing references to these objects
    # is added as an attribute to the std::map.
    arg_map = ROOT.std.map[get_template_args(arg_dict)]()
    arg_map.keepalive = list()
    for key, value in arg_dict.items():
        arg_map.keepalive.append(value)
        arg_map.emplace(key, value)

    return arg_map


def _decaytype_string_to_enum(caller, kwargs):
    """Helper function to pythonize DecayType enums and check for enum value names."""
    type_key = "type"

    if type_key in kwargs:
        val = kwargs[type_key]
        if isinstance(val, str):
            try:
                kwargs[type_key] = getattr(caller.__class__, val)
            except AttributeError as error:
                raise ValueError(
                    "Unsupported decay type passed to "
                    + caller.__class__.__name__
                    + ". Supported decay types are : 'SingleSided', 'DoubleSided', 'Flipped'"
                )
            except Exception as exception:
                raise exception

    return kwargs
