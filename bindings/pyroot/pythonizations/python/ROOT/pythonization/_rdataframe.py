# Author: Stefan Wunsch, Massimiliano Galli CERN  02/2019

################################################################################
# Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization

# functools.partial does not add the self argument
# this is done by functools.partialmethod which is
# introduced only in Python 3.4
try:
    from functools import partialmethod
except ImportError:
    from functools import partial

    class partialmethod(partial):
        def __get__(self, instance, owner):
            if instance is None:
                return self
            return partial(self.func, instance, *(self.args or ()), **(self.keywords or {}))


def RDataFrameAsNumpy(df, columns=None, exclude=None):
    """Read-out the RDataFrame as a collection of numpy arrays.

    The values of the dataframe are read out as numpy array of the respective type
    if the type is a fundamental type such as float or int. If the type of the column
    is a complex type, such as your custom class or a std::array, the returned numpy
    array contains Python objects of this type interpreted via PyROOT.

    Be aware that reading out custom types is much less performant than reading out
    fundamental types, such as int or float, which are supported directly by numpy.

    The reading is performed in multiple threads if the implicit multi-threading of
    ROOT is enabled.

    Note that this is an instant action of the RDataFrame graph and will trigger the
    event-loop.

    Parameters:
        columns: If None return all branches as columns, otherwise specify names in iterable.
        exclude: Exclude branches from selection.

    Returns:
        dict: Dict with column names as keys and 1D numpy arrays with content as values
    """
    # Sanitize input arguments
    if isinstance(columns, str):
        raise TypeError("The columns argument requires a list of strings")
    if isinstance(exclude, str):
        raise TypeError("The exclude argument requires a list of strings")

    # Import numpy and numpy.array derived class lazily
    try:
        import numpy
        from ROOT.pythonization._rdf_utils import ndarray
    except:
        raise ImportError("Failed to import numpy during call of RDataFrame.AsNumpy.")

    # Find all column names in the dataframe if no column are specified
    if not columns:
        columns = [str(c) for c in df.GetColumnNames()]

    # Exclude the specified columns
    if exclude == None:
        exclude = []
    columns = [col for col in columns if not col in exclude]

    # Register Take action for each column
    result_ptrs = {}
    for column in columns:
        column_type = df.GetColumnType(column)
        result_ptrs[column] = df.Take[column_type](column)

    # Convert the C++ vectors to numpy arrays
    py_arrays = {}
    for column in columns:
        cpp_reference = result_ptrs[column].GetValue()
        if hasattr(cpp_reference, "__array_interface__"):
            tmp = numpy.asarray(cpp_reference) # This adopts the memory of the C++ object.
            py_arrays[column] = ndarray(tmp, result_ptrs[column])
        else:
            tmp = numpy.empty(len(cpp_reference), dtype=numpy.object)
            for i, x in enumerate(cpp_reference):
                tmp[i] = x # This creates only the wrapping of the objects and does not copy.
            py_arrays[column] = ndarray(tmp, result_ptrs[column])

    return py_arrays


def _histo_profile(self, fixed_args, *args):
    # Check wheter the user called one of the HistoXD or ProfileXD methods
    # of RDataFrame with a tuple as first argument; in that case,
    # extract the tuple items to construct a model object and call the
    # original implementation of the method with that object.

    # Parameters:
    # self: instantiation of RDataFrame
    # fixed_args: tuple containing the original name of the method being
    # pythonised and the class of the model object to construct
    # args: arguments passed by the user when he calls e.g Histo1D

    original_method_name, model_class = fixed_args

    # Get the "original" method of the RDataFrame instantiation
    original_method = getattr(self, original_method_name)

    if args and isinstance(args[0], tuple):
        # Construct the model with the elements of the tuple
        # as arguments
        model = model_class(*args[0])
        # Call the original implementation of the method
        # with the model as first argument
        if len(args) > 1:
            res = original_method(model, *args[1:])
        else:
            # Covers the case of the overloads with only model passed
            # as argument
            res = original_method(model)
    # If the first argument is not a tuple, nothing to do, just call
    # the original implementation
    else:
        res = original_method(*args)

    return res


@pythonization()
def pythonize_rdataframe(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class

    if name.startswith("ROOT::RDataFrame<") or name.startswith("ROOT::RDF::RInterface<"):
        from cppyy.gbl.ROOT import RDF

        # Add asNumpy feature
        klass.AsNumpy = RDataFrameAsNumpy

        # Replace the implementation of the following RDF methods
        # to convert a tuple argument into a model object
        methods_with_TModel = {
                'Histo1D' : RDF.TH1DModel,
                'Histo2D' : RDF.TH2DModel,
                'Histo3D' : RDF.TH3DModel,
                'Profile1D' : RDF.TProfile1DModel,
                'Profile2D' : RDF.TProfile2DModel
                }

        # Do e.g.:
        # klass._OriginalHisto1D = klass.Histo1D
        # klass.Histo1D = TH1DModel
        for method_name, model_class in methods_with_TModel.items():
            original_method_name = '_Original' + method_name
            setattr(klass, original_method_name, getattr(klass, method_name))
            # Fixed arguments to construct a partialmethod
            fixed_args = (original_method_name, model_class)
            # Replace the original implementation of the method
            # by a generic function _histo_profile with
            # (original_method_name, model_class) as fixed argument
            setattr(klass, method_name, partialmethod(_histo_profile, fixed_args))

    return True
