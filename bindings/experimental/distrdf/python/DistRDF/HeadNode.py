import glob
import logging
import warnings

import ROOT

from DistRDF import Node
from DistRDF import Ranges

logger = logging.getLogger(__name__)


def get_headnode(npartitions, *args):
    """
    A factory for different kinds of head nodes of the RDataFrame computation
    graph, depending on the arguments to the RDataFrame constructor. Currently
    can return a TreeHeadNode or an EmptySourceHeadNode. Parses the arguments and
    compares them against the possible RDataFrame constructors.
    """

    # Early check that arguments are accepted by RDataFrame
    try:
        ROOT.RDataFrame(*args)
    except TypeError:
        raise TypeError(("The arguments provided are not accepted by any RDataFrame constructor. "
                         "See the RDataFrame documentation for the accepted constructor argument types."))

    firstarg = args[0]
    if isinstance(firstarg, int):
        # RDataFrame(ULong64_t numEntries)
        return EmptySourceHeadNode(npartitions, firstarg)
    elif isinstance(firstarg, (ROOT.TTree, str)):
        # RDataFrame(std::string_view treeName, filenameglob, defaultBranches = {})
        # RDataFrame(std::string_view treename, filenames, defaultBranches = {})
        # RDataFrame(std::string_view treeName, dirPtr, defaultBranches = {})
        # RDataFrame(TTree &tree, const ColumnNames_t &defaultBranches = {})
        return TreeHeadNode(npartitions, *args)
    else:
        raise RuntimeError(
            ("First argument {} of type {} is not recognised as a supported "
                "argument for distributed RDataFrame. Currently only TTree/Tchain "
                "based datasets or datasets created from a number of entries "
                "can be processed distributedly.").format(firstarg, type(firstarg)))


class EmptySourceHeadNode(Node.Node):
    """
    The head node of a computation graph where the RDataFrame data source is
    empty and a number of sequential entries will be created at runtime. This
    head node is responsible for the following RDataFrame constructor::

        RDataFrame(ULong64_t numEntries)

    Attributes:
        nentries (int): The number of sequential entries the RDataFrame will create.

        npartitions (int): The number of partitions the dataset will be split in
            for distributed execution.
    """

    def __init__(self, npartitions, nentries):
        """
        Creates a new RDataFrame instance for the given arguments.

        Args:
            nentries (int): The number of entries this RDataFrame will process.

            npartitions (int): The number of partitions the dataset will be
                split in for distributed execution.
        """
        super(EmptySourceHeadNode, self).__init__(None, None)

        self.nentries = nentries
        self.npartitions = npartitions

    def build_ranges(self):
        """Build the ranges for this dataset."""
        # Empty datasets cannot be processed distributedly
        if not self.nentries:
            raise RuntimeError(
                ("Cannot build a distributed RDataFrame with zero entries. "
                 "Distributed computation will fail. "))
        # TODO: This shouldn't be triggered if entries == 1. The current minimum
        # amount of partitions is 2. We need a robust reducer that smartly
        # becomes no-op if npartitions == 1 to avoid this.
        if self.npartitions > self.nentries:
            # Restrict 'npartitions' if it's greater than 'nentries'
            msg = ("Number of partitions {0} is greater than number of entries {1} "
                   "in the dataframe. Using {1} partition(s)".format(self.npartitions, self.nentries))
            warnings.warn(msg, UserWarning, stacklevel=2)
            self.npartitions = self.nentries
        return Ranges.get_balanced_ranges(self.nentries, self.npartitions)


class TreeHeadNode(Node.Node):
    """
    The head node of a computation graph where the RDataFrame data source is
    a TTree. This head node is responsible for the following RDataFrame constructor::

        RDataFrame(std::string_view treeName, std::string_view filenameglob, const ColumnNames_t &defaultBranches = {})
        RDataFrame(std::string_view treename, const std::vector<std::string> &fileglobs, const ColumnNames_t &defaultBranches = {})
        RDataFrame(std::string_view treeName, TDirectory *dirPtr, const ColumnNames_t &defaultBranches = {})
        RDataFrame(TTree &tree, const ColumnNames_t &defaultBranches = {})

    Attributes:
        args (iterable): An iterable of arguments that were provided to construct
            the RDataFrame object.

        npartitions (int): The number of partitions the dataset will be split in
            for distributed execution.
    """

    def __init__(self, npartitions, *args):
        """
        Creates a new RDataFrame instance for the given arguments.

        Args:
            *args (iterable): Iterable with the arguments to the RDataFrame constructor.

            npartitions (int): Keyword argument with the number of partitions
                the dataset will be split in for distributed execution.
        """
        super(TreeHeadNode, self).__init__(None, None)

        # Keep the arguments to parse them when needed
        # TODO: Instead of storing the args here and parsing them various times
        # in different class methods, parse them only once in this function and
        # directly store the appropriate data members
        self.args = args
        self.npartitions = npartitions

    # TODO: Decide whether to remove/change the property or the getter
    @property
    def tree(self):
        """Tree instance if present."""
        return self.get_tree()

    # TODO: Decide whether to remove/change the property or the getter
    @property
    def treename(self):
        """Name of the tree."""
        return self.get_treename()

    # TODO: Decide whether to remove/change the property or the getter
    @property
    def nentries(self):
        """Entries of this dataset."""
        return self.get_num_entries()

    # TODO: Decide whether to remove/change the property or the getter
    @property
    def defaultbranches(self):
        """Default branches selected by the user in the constructor."""
        return self.get_branches()

    # TODO: Decide whether to remove/change the property or the getter
    @property
    def inputfiles(self):
        """List of input files of the dataset."""
        return self.get_inputfiles()

    # TODO: Decide whether to remove/change the property or the getter
    @property
    def friendinfo(self):
        """Information about friend trees of the dataset."""
        return self._get_friend_info()

    def get_branches(self):
        """Gets list of default branches if passed by the user."""
        # ROOT Constructor:
        # RDataFrame(TTree& tree, defaultBranches = {})
        if len(self.args) == 2 and isinstance(self.args[0], ROOT.TTree):
            return self.args[1]
        # ROOT Constructors:
        # RDataFrame(treeName, filenameglob, defaultBranches = {})
        # RDataFrame(treename, filenames, defaultBranches = {})
        # RDataFrame(treeName, dirPtr, defaultBranches = {})
        if len(self.args) == 3:
            return self.args[2]

        return None

    # TODO: This function is very costly! Make sure to use only once at most, 
    # or think if it can be avoided altogether.
    def get_num_entries(self):
        """
        Gets the number of entries in the given dataset.

        Returns:
            int: This is the computed number of entries in the input dataset.

        """
        # First argument can be a string or a TTree/TChain
        firstarg = self.args[0]
        if isinstance(firstarg, ROOT.TTree):
            # If the argument is a TTree or TChain,
            # get the number of entries from it.
            return firstarg.GetEntries()
        elif isinstance(firstarg, str):
            # Construct a ROOT.TChain object with the name of the tree
            chain = ROOT.TChain(firstarg)

            # Add the list of input files to the chain
            for fname in self.inputfiles:
                chain.Add(fname)

            return chain.GetEntries()

    def get_treename(self):
        """
        Get name of the input tree. If the user supplied a string as first
        argument, that corresponds to the tree name. If the first argument is a
        TTree or TChain, it retrieves the first full path from the paths
        returned by ROOT::Internal::TreeUtils::GetTreeFullPaths function.

        Returns:
            (str): Name of the tree. 

        """
        firstarg = self.args[0]
        if isinstance(firstarg, ROOT.TTree):
            treefullpaths = ROOT.Internal.TreeUtils.GetTreeFullPaths(firstarg)
            return treefullpaths[0]
        elif isinstance(firstarg, str):
            # First argument was the name of the tree
            return firstarg

    def get_tree(self):
        """
        Get ROOT.TTree instance given as constructor argument.

        Returns:
            (ROOT.TTree, None): instance of the tree used to instantiate the
            RDataFrame, or `None` if another object was used. ROOT.Tchain
            inherits from ROOT.TTree so that can be the return value as well.
        """
        first_arg = self.args[0]
        if isinstance(first_arg, ROOT.TTree):
            return first_arg

        return None

    def get_inputfiles(self):
        """
        Get list of input files of the dataset. If the user supplied a TTree or
        TChain, returns the result of calling the
        ROOT::Internal::TreeUtils::GetFileNamesFromTree function on that tree.
        If one of the other tree based constructor arguments of RDataFrame was
        supplied, this function will parse them and return the name(s) of the
        corresponding file(s).

        Returns:
            (list[str]): List of input files of the dataset.

        """
        if isinstance(self.args[0], ROOT.TTree):
            # Extract file names from a given TTree or TChain
            filenames = ROOT.Internal.TreeUtils.GetFileNamesFromTree(self.args[0])
            return [str(filename) for filename in filenames]

        if len(self.args) > 1:
            # Second argument can be:
            # 1. A simple string representing a file or a glob of files
            # 2. A vector of filename strings
            # 3. A pointer to a TDirectory (i.e. the file of the dataset)
            secondarg = self.args[1]
            if isinstance(secondarg, str):
                # Expand globbing excluding remote files
                remote_prefixes = ("root:", "http:", "https:")
                if not secondarg.startswith(remote_prefixes):
                    return glob.glob(secondarg)
                else:
                    return [secondarg]
            elif isinstance(secondarg, (list, ROOT.std.vector("string"))):
                # Make sure this returns a list of Python strings
                return [str(filename) for filename in secondarg]
            elif isinstance(secondarg, ROOT.TDirectory):
                # Return the name of the file
                return [str(secondarg.GetName())]

    def _get_friend_info(self):
        """
        Retrieves information about friend trees of the input tree.

        Returns:
            (ROOT.Internal.TreeUtils.RFriendInfo): A C++ struct holding
                information about friend trees. If the user did not supply a
                TTree or TChain as input to the constructor, returns None.
        """
        if isinstance(self.args[0], ROOT.TTree):
            return ROOT.Internal.TreeUtils.GetFriendInfo(self.args[0])
        else:
            return None

    def build_ranges(self):
        """Build the ranges for this dataset."""
        # Empty datasets cannot be processed distributedly
        if not self.nentries:
            raise RuntimeError(
                ("Cannot build a distributed RDataFrame with zero entries. "
                 "Distributed computation will fail. "))

        logger.debug("Building ranges for tree %s with the "
                     "following input files:\n%s", self.treename, self.inputfiles)

        # Retrieve a tuple of clusters for all files of the tree
        clustersinfiles = Ranges.get_clusters(self.treename, self.inputfiles)
        numclusters = len(clustersinfiles)

        # TODO: This shouldn't be triggered if len(clustersinfiles) == 1. The
        # current minimum amount of partitions is 2. We need a robust reducer
        # that smartly becomes no-op if npartitions == 1 to avoid this.
        # Restrict `npartitions` if it's greater than clusters of the dataset
        if self.npartitions > numclusters:
            msg = ("Number of partitions is greater than number of clusters "
                   "in the dataset. Using {} partition(s)".format(numclusters))
            warnings.warn(msg, UserWarning, stacklevel=2)
            self.npartitions = numclusters

        logger.debug("%s clusters will be split along %s partitions.",
                     numclusters, self.npartitions)
        return Ranges.get_clustered_ranges(clustersinfiles, self.npartitions, self.friendinfo)
