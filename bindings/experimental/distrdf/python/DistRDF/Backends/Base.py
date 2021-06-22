# @author Vincenzo Eduardo Padulano
#  @author Enric Tejedor
#  @date 2021-02

################################################################################
# Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

import functools
from abc import ABCMeta
from abc import abstractmethod

import ROOT
from DistRDF.Backends import Utils

# Abstract class declaration
# This ensures compatibility between Python 2 and 3 versions, since in
# Python 2 there is no ABC class
ABC = ABCMeta("ABC", (object,), {})


class BaseBackend(ABC):
    """
    Base class for RDataFrame distributed backends.

    Attributes:
        supported_operations (list): List of operations supported by the
            backend.
        initialization (function): Store user's initialization method, if
            defined.
        headers (list): List of headers that need to be declared for the
            analysis.
        shared_libraries (list): List of shared libraries needed for the
            analysis.
    """

    supported_operations = [
        "AsNumpy",
        "Count",
        "Define",
        "Fill",
        "Filter",
        "Graph",
        "Histo1D",
        "Histo2D",
        "Histo3D",
        "Max",
        "Mean",
        "Min",
        "Profile1D",
        "Profile2D",
        "Profile3D",
        "Snapshot",
        "Sum"
    ]

    initialization = staticmethod(lambda: None)

    headers = set()
    shared_libraries = set()

    @classmethod
    def register_initialization(cls, fun, *args, **kwargs):
        """
        Convert the initialization function and its arguments into a callable
        without arguments. This callable is saved on the backend parent class.
        Therefore, changes on the runtime backend do not require users to set
        the initialization function again.

        Args:
            fun (function): Function to be executed.

            *args (list): Variable length argument list used to execute the
                function.

            **kwargs (dict): Keyword arguments used to execute the function.
        """
        cls.initialization = functools.partial(fun, *args, **kwargs)
        fun(*args, **kwargs)

    def check_supported(self, operation_name):
        """
        Checks if a given operation is supported
        by the given backend.

        Args:
            operation_name (str): Name of the operation to be checked.

        Raises:
            Exception: This happens when `operation_name` doesn't exist
            the `supported_operations` instance attribute.
        """
        if operation_name not in self.supported_operations:
            raise Exception(
                "The current backend doesn't support \"{}\" !"
                .format(operation_name)
            )

    def execute(self, generator):
        """
        Executes an RDataFrame computation graph on a distributed backend.

        Args:
            generator (ComputationGraphGenerator): A factory object for a
                computation graph. Its ``get_callable`` method will return a
                function responsible for creating the computation graph of a
                given RDataFrame object and a range of entries. The range is
                needed for the `Snapshot` operation.
        """
        headnode = generator.headnode
        computation_graph_callable = generator.get_callable()

        # Avoid having references to the instance inside the mapper
        initialization = self.initialization

        # Build the ranges for the current dataset
        ranges = headnode.build_ranges()
        headnode_type = headnode.TYPE

        def mapper(current_range):
            """
            Triggers the event-loop and executes all
            nodes in the computational graph using the
            callable.

            Args:
                current_range (Range): A Range named tuple, representing the
                    range of entries to be processed, their input files and
                    information about friend trees.

            Returns:
                list: This respresents the list of (mergeable)values of all
                action nodes in the computational graph.
            """
            import ROOT

            print(current_range)
            # We have to decide whether to do this in Dist or in subclasses
            # Utils.declare_headers(worker_includes)  # Declare headers if any
            # Run initialization method to prepare the worker runtime
            # environment
            initialization()

            # Build rdf
            start = int(current_range.start)
            end = int(current_range.end)

            # Check the type of the headnode. Each headnode will create Range
            # objects of different kinds and in the distributed mapper will be
            # needed to create the RDataFrame object according to the constructor
            # arguments. The following two approaches cannot work
            # 1. Having separate helper functions to create the RDataFrames,
            #    for example `rdf = headnode.get_rdf_for_mapper()`, because
            #    RDataFrame is not serializable (it misses a default constructor)
            # 2. Checking for the type of the Range object, e.g.
            #    `if isinstance(current_range, Ranges.TreeRange`, because the
            #    type of a namedtuple is actually `collections.TreeRange` but
            #    it's not sufficient to `import collections` at the beginning
            #    of this file.
            # For now the easiest and most lightweight approach is to store a
            # string with a certain name unique to the headnode and use that
            # as a discriminant.
            if headnode_type == "TREE":
                # Build TChain of files for this range:
                chain = ROOT.TChain(current_range.treename)
                for f in current_range.treefilenames:
                    chain.Add(str(f))

                # We assume 'end' is exclusive
                chain.SetCacheEntryRange(start, end)

                # Gather information about friend trees
                if current_range.friendnamesalias and current_range.friendfilenames and current_range.friendchainsubnames:
                    # Zip together the information about friend trees. Each element if of the form:
                    # (name, alias), (file1.root, file2.root, ...), (subname1, subname2, ...)
                    tree_files_names = zip(
                        current_range.friendnamesalias,
                        current_range.friendfilenames,
                        current_range.friendchainsubnames
                    )
                    for (friend_name, friend_alias), friend_filenames, friendchainsubnames in tree_files_names:
                        # Start a TChain with the current friend treename
                        friend_chain = ROOT.TChain(friend_name)
                        # Add each corresponding file to the TChain
                        for filename, chainsubname in zip(friend_filenames, friendchainsubnames):
                            # filename is an std::string, Add accepts a const char*
                            if not chainsubname:
                                friend_chain.Add(filename)
                            else:
                                friend_chain.Add(filename+"/"+chainsubname)

                        # Set cache on the same range as the parent TChain
                        friend_chain.SetCacheEntryRange(start, end)
                        # Finally add friend TChain to the parent (with alias)
                        chain.AddFriend(friend_chain, friend_alias)

                if current_range.defaultbranches is not None:
                    rdf = ROOT.RDataFrame(chain, current_range.defaultbranches)
                else:
                    rdf = ROOT.RDataFrame(chain)
            elif headnode_type == "ENTRIES":
                rdf = ROOT.RDataFrame(current_range.nentries)

            # Output of the callable
            resultptr_list = computation_graph_callable(
                rdf, rdf_range=current_range)

            mergeables = [
                resultptr  # Here resultptr is already the result value
                if isinstance(resultptr, (dict, list))
                else ROOT.ROOT.Detail.RDF.GetMergeableValue(resultptr)
                for resultptr in resultptr_list
            ]
            return mergeables

        def reducer(mergeables_out, mergeables_in):
            """
            Merges two given lists of values that were
            returned by the mapper function for two different
            ranges.

            Args:
                mergeables_out (list): A list of computed (mergeable)values for
                    a given entry range in a dataset. The elements of this list
                    will be updated with the information contained in the
                    elements of the other argument list.

                mergeables_in (list): A list of computed (mergeable)values for
                    a given entry range in a dataset.

            Returns:
                list: The list of updated (mergeable)values.
            """

            import ROOT

            # We still need the list index to modify results of `Snapshot` and
            # `AsNumpy` in place.
            for index, (mergeable_out, mergeable_in) in enumerate(
                    zip(mergeables_out, mergeables_in)):
                # Create a global list with all the files of the partial
                # snapshots.
                if isinstance(mergeable_out, list):
                    mergeables_out[index].extend(mergeable_in)

                # Concatenate the partial numpy arrays along the same key of
                # the dictionary.
                elif isinstance(mergeable_out, dict):
                    # Import numpy lazily
                    try:
                        import numpy
                    except ImportError:
                        raise ImportError("Failed to import numpy during distributed RDataFrame reduce step.")
                    mergeables_out[index] = {
                        key: numpy.concatenate([mergeable_out[key],
                                                mergeable_in[key]])
                        for key in mergeable_out
                    }

                # The `MergeValues` function modifies the arguments in place
                # so there's no need to access the list elements.
                else:
                    ROOT.ROOT.Detail.RDF.MergeValues(
                        mergeable_out, mergeable_in)

            return mergeables_out

        # Values produced after Map-Reduce
        values = self.ProcessAndMerge(ranges, mapper, reducer)
        # List of action nodes in the same order as values
        nodes = generator.get_action_nodes()

        # Set the value of every action node
        for node, value in zip(nodes, values):
            if node.operation.name == "Snapshot":
                # Retrieve treename from operation args and start TChain
                snapshot_treename = node.operation.args[0]
                snapshot_chain = ROOT.TChain(snapshot_treename)
                # Add partial snapshot files to the chain
                for filename in value:
                    snapshot_chain.Add(filename)
                # Create a new rdf with the chain and return that to user
                node.value = self.make_dataframe(snapshot_chain)
            elif node.operation.name == "AsNumpy":
                node.value = value
            else:
                node.value = value.GetValue()

    @abstractmethod
    def ProcessAndMerge(self, ranges, mapper, reducer):
        """
        Subclasses must define how to run map-reduce functions on a given
        backend.
        """
        pass

    @abstractmethod
    def distribute_unique_paths(self, paths):
        """
        Subclasses must define how to send all files needed for the analysis
        (like headers and libraries) to the workers.
        """
        pass

    def optimize_npartitions(self, npartitions):
        """
        Distributed backends may optimize the number of partitions of the
        current dataset or leave it as it is.
        """
        return npartitions

    def distribute_files(self, files_paths):
        """
        Sends to the workers the generic files needed by the user.

        Args:
            files_paths (str, iter): Paths to the files to be sent to the
                distributed workers.
        """
        files_to_distribute = set()

        if isinstance(files_paths, str):
            files_to_distribute.update(
                Utils.get_paths_set_from_string(files_paths))
        else:
            for path_string in files_paths:
                files_to_distribute.update(
                    Utils.get_paths_set_from_string(path_string))

        self.distribute_unique_paths(files_to_distribute)

    def distribute_headers(self, headers_paths):
        """
        Includes the C++ headers to be declared before execution.

        Args:
            headers_paths (str, iter): A string or an iterable (such as a
                list, set...) containing the paths to all necessary C++ headers
                as strings. This function accepts both paths to the headers
                themselves and paths to directories containing the headers.
        """
        headers_to_distribute = set()

        if isinstance(headers_paths, str):
            headers_to_distribute.update(
                Utils.get_paths_set_from_string(headers_paths))
        else:
            for path_string in headers_paths:
                headers_to_distribute.update(
                    Utils.get_paths_set_from_string(path_string))

        # Distribute header files to the workers
        self.distribute_unique_paths(headers_to_distribute)

        # Declare headers locally
        Utils.declare_headers(headers_to_distribute)

        # Finally, add everything to the includes set
        self.headers.update(headers_to_distribute)

    def distribute_shared_libraries(self, shared_libraries_paths):
        """
        Includes the C++ shared libraries to be declared before execution. If
        any pcm file is present in the same folder as the shared libraries, the
        function will try to retrieve them and distribute them.

        Args:
            shared_libraries_paths (str, iter): A string or an iterable (such as
                a list, set...) containing the paths to all necessary C++ shared
                libraries as strings. This function accepts both paths to the
                libraries themselves and paths to directories containing the
                libraries.
        """
        libraries_to_distribute = set()
        pcm_to_distribute = set()

        if isinstance(shared_libraries_paths, str):
            pcm_to_distribute, libraries_to_distribute = (
                Utils.check_pcm_in_library_path(shared_libraries_paths))
        else:
            for path_string in shared_libraries_paths:
                pcm, libraries = Utils.check_pcm_in_library_path(
                    path_string
                )
                libraries_to_distribute.update(libraries)
                pcm_to_distribute.update(pcm)

        # Distribute shared libraries and pcm files to the workers
        self.distribute_unique_paths(libraries_to_distribute)
        self.distribute_unique_paths(pcm_to_distribute)

        # Include shared libraries locally
        Utils.declare_shared_libraries(libraries_to_distribute)

        # Finally, add everything to the includes set
        self.shared_libraries.update(libraries_to_distribute)

    @abstractmethod
    def make_dataframe(self, *args, **kwargs):
        """
        Distributed backends have to take care of creating an RDataFrame object
        that can run distributedly.
        """
