import glob
import logging
import warnings

import ROOT

from DistRDF import Node
from DistRDF import Ranges

logger = logging.getLogger(__name__)


class Factory(object):
    """
    A factory for different kinds of root nodes of the RDataFrame computation
    graph, depending on the arguments to the RDataFrame constructor. Currently
    can create a TreeHeadNode or an EntriesHeadNode.
    """

    @staticmethod
    def get_headnode(*args):
        """
        Parses the arguments and compares them against the possible RDataFrame
        constructors
        """

        # Early check that arguments are accepted by RDataFrame
        ROOT.RDataFrame(*args)

        firstarg = args[0]
        if isinstance(firstarg, int):
            # RDataFrame(ULong64_t numEntries)
            return EntriesHeadNode(firstarg)
        elif isinstance(firstarg, (ROOT.TTree, str)):
            # RDataFrame(std::string_view treeName, filenameglob, defaultBranches = {})
            # RDataFrame(std::string_view treename, filenames, defaultBranches = {})
            # RDataFrame(std::string_view treeName, dirPtr, defaultBranches = {})
            # RDataFrame(TTree &tree, const ColumnNames_t &defaultBranches = {})
            return TreeHeadNode(*args)
        else:
            raise RuntimeError(
                ("First argument {} of type {} is not recognised as a supported "
                 "argument for distributed RDataFrame. Currently only TTree/Tchain "
                 "based datasets or datasets created from a number of entries "
                 "can be processed distributedly.").format(firstarg, type(firstarg)))


class EntriesHeadNode(Node.Node):
    """
    The head node of a computation graph where the RDataFrame is created from
    a number of entries.
    """

    TYPE = "ENTRIES"

    def __init__(self, nentries):
        """
        Creates a new RDataFrame instance for the given arguments.

        Args:
            nentries (int): The number of entries this RDataFrame will process.
        """
        super(EntriesHeadNode, self).__init__(None, None)

        self.nentries = nentries
        # Set at creation of the dataframe, might be optimized by the backend
        # in optimize_partitions
        self.npartitions = 2

    def build_ranges(self):
        """Build the ranges for this dataset."""
        # Empty datasets cannot be processed distributedly
        if not self.nentries:
            raise RuntimeError(
                ("Cannot build a distributed RDataFrame with zero entries. "
                 "Distributed computation will fail. "))
        if self.npartitions > self.nentries:
            # Restrict 'npartitions' if it's greater than 'nentries'
            msg = ("Number of partitions {0} is greater than number of entries {1} "
                   "in the dataframe. Using {1} partition(s)".format(self.npartitions, self.nentries))
            warnings.warn(msg, UserWarning, stacklevel=2)
            self.npartitions = self.nentries
        return Ranges.get_balanced_ranges(self.nentries, self.npartitions)


class TreeHeadNode(Node.Node):
    """
    The Python equivalent of ROOT C++'s
    RDataFrame class.

    Attributes:
        args (list): A list of arguments that were provided to construct
            the RDataFrame object.


    DistRDF's RDataFrame constructor accepts the same arguments as the ROOT's
    RDataFrame constructor (see
    `RDataFrame <https://root.cern/doc/master/classROOT_1_1RDataFrame.html>`_)
    """

    TYPE = "TREE"

    def __init__(self, *args):
        """
        Creates a new RDataFrame instance for the given arguments.

        Args:
            *args (list): Variable length argument list to construct the
                RDataFrame object.
        """
        super(TreeHeadNode, self).__init__(None, None)

        # Keep the arguments to parse them when needed
        self.args = args
        # Set at creation of the dataframe, might be optimized by the backend
        # in optimize_partitions
        self.npartitions = 2

    @property
    def tree(self):
        """Tree instance if present"""
        return self.get_tree()

    @property
    def treename(self):
        """Name of the tree"""
        return self.get_treename()

    @property
    def nentries(self):
        """Entries of this dataset"""
        return self.get_num_entries()

    @property
    def defaultbranches(self):
        """Default branches selected by the user in the constructor"""
        return self.get_branches()

    @property
    def inputfiles(self):
        """List of input files of the dataset if any"""
        return self.get_inputfiles()

    @property
    def friendinfo(self):
        """Retrieve information about friend trees"""
        return self.get_friendinfo()

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
        else:
            raise RuntimeError("Could not retrieve number of entries from the dataset.")

    def get_treename(self):
        """
        Get name of the TTree.

        Returns:
            (str, None): Name of the TTree, or :obj:`None` if there is no tree.

        """
        firstarg = self.args[0]
        if isinstance(firstarg, ROOT.TTree):
            treefullpaths = ROOT.Internal.TreeUtils.GetTreeFullPaths(firstarg)
            return treefullpaths[0]
        elif isinstance(firstarg, str):
            # First argument was the name of the tree
            return firstarg
        else:
            raise RuntimeError("Could not find name of the input tree.")

    def get_tree(self):
        """
        Get ROOT.TTree instance used as an argument to DistRDF.RDataFrame()

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
        Get list of input files.

        This list can be extracted from a given TChain or from the list of
        arguments.

        Returns:
            (list, None): List of files, or None if there are no input files.

        """
        firstarg = self.args[0]
        if isinstance(firstarg, ROOT.TTree):
            # Extract file names from a given TChain
            filenames = ROOT.Internal.TreeUtils.GetFileNamesFromTree(firstarg)
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

        # We should be able to retrieve file names from the input arguments.
        # Otherwise, error out since distributed processing of a tree with no
        # input files is not supported.
        raise RuntimeError("Could not find input files of the tree.")

    def get_friendinfo(self):
        """
        Retrieve information about the friends of the tree used to construct
        this RDataFrame
        Returns:
            (tuple): The list of names and aliases of the friends and their
                corresponding filenames. If the argument to the constructor is
                not a TTree/TChain the elements of the tuple are None.
        """
        # Reconstruct friendinfo from the input tree
        firstarg = self.args[0]
        if isinstance(firstarg, ROOT.TTree):
            # RDataFrame(TTree &tree, const ColumnNames_t &defaultBranches = {})
            # The first argument to the constructor is a TTree or TChain
            treefriendinfo = ROOT.Internal.TreeUtils.GetFriendInfo(firstarg)

            # Convert to tuples to be able to hash them and cache the clustered ranges
            friendnamesalias = tuple((str(namealias.first), str(namealias.second))
                                     for namealias in treefriendinfo.fFriendNames)
            friendfilenames = tuple(tuple(str(filename) for filename in filenames)
                                    for filenames in treefriendinfo.fFriendFileNames)
            friendchainsubnames = tuple(tuple(str(subname) for subname in chainsubnames)
                                    for chainsubnames in treefriendinfo.fFriendChainSubNames)
            return friendnamesalias, friendfilenames, friendchainsubnames
        else:
            return None, None, None

    def build_ranges(self):
        """Build the ranges for this dataset."""

        # Empty trees cannot be processed distributedly
        if not self.nentries:
            raise RuntimeError(
                ("No entries in the TTree. "
                 "Distributed computation will fail. "
                 "Please make sure your dataset is not empty."))

        treename = self.treename
        inputfiles = self.inputfiles
        friendnamesalias, friendfilenames, friendchainsubnames = self.friendinfo
        defaultbranches = self.defaultbranches

        logger.debug("Building ranges for tree %s with the "
                     "following input files:\n%s", treename, inputfiles)

        # Retrieve a tuple of clusters for all files of the tree
        clustersinfiles = Ranges.get_clusters(treename, tuple(inputfiles))
        numclusters = len(clustersinfiles)

        # Restrict `npartitions` if it's greater than clusters of the dataset
        if self.npartitions > numclusters:
            msg = ("Number of partitions is greater than number of clusters "
                   "in the dataset. Using {} partition(s)".format(numclusters))
            warnings.warn(msg, UserWarning, stacklevel=2)
            self.npartitions = numclusters

        logger.debug("%s clusters will be split along %s partitions.",
                     numclusters, self.npartitions)
        return Ranges.get_clustered_ranges(tuple(clustersinfiles), self.npartitions, treename, friendnamesalias,
                                           friendfilenames, friendchainsubnames, defaultbranches)
